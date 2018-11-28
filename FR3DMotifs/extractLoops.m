%==========================================================================
% Extracts all RNA 3D motifs of a specified loopType from a set of PDB files.
% Input: pdb id, motif loopType (IL, HL, J3, J4 etc).
% Returns loops as a structure array or 0 if no loops were found.
% All motifs, including noncontiguous and ncWW are stored. Postprocessing
% should be done separately.
%==========================================================================
function [Loops, L, err_msg] = extractLoops(pdb_id,loopType)

    global loop_type;

    try
        loop_type = upper(loopType);
        Loops     = 0;
        L         = 0;
        err_msg   = '';

        if strcmp(loopType,'HL')
            Loops = aFlankHairpinSearch(pdb_id);
        elseif strcmp(loopType,'IL')
            Loops = aCombinedInternalLoopSearch(pdb_id);
        elseif strfind(loopType,'J')
            Loops = aFlankJunctionSearch(pdb_id,loopType);
        end

        L = length(Loops);
        aCleanUp;

    catch err
        err_msg = sprintf('Error "%s" in extractLoops on line %i (%s)\n', ...
                           err.message, err.stack.line, pdb_id);
        disp(err_msg);
    end

end

%==========================================================================
% Create FR3D `File` structures for each loop.
%==========================================================================
function [Loops] = aSearchToLoops(Search)

    global loop_type;

    for i = 1:length(Search.Candidates(:,1))

        pdbfile = Search.Candidates(i,end);
        File    = Search.File(pdbfile);
        [range,breaks,ranges,nwcranges] = aGetRange(Search.Candidates(i,:));
        File.chain_breaks = breaks;

        fn = fieldnames(File);
        for j = 1:length(fn)
            [r,c] = size(File.(fn{j}));
            if r == c && r == File.NumNT
                File.(fn{j}) = File.(fn{j})(range,range);
            end
        end

        File.NT    = File.NT(range);
        File.NumNT = length(range);

        % to save space
        if isfield(File,'Het'), File.Het = []; end
        if isfield(File,'AA'),  File.AA  = []; end

        N = length(File.NT);
        ids = cell(1,N);
        for j = 1:N
            ids{j} = [File.NT(j).ID ','];
        end
        ids = sort(ids);
        File.AllLoops_table.full_id = [ids{:}];
        File.AllLoops_table.full_id = File.AllLoops_table.full_id(1:end-1);

        File.AllLoops_table.loop_name = aGetLoopName(File);

        nt_seq   = [File.NT.Base];
        if strcmp(loop_type,'HL')
            File.AllLoops_table.seq   = nt_seq;
            File.AllLoops_table.r_seq = '';
            File.AllLoops_table.nwc   = File.AllLoops_table.seq(2:end-1);
            File.AllLoops_table.r_nwc = '';
            File.AllLoops_table.sequence{1}    = nt_seq;
            File.AllLoops_table.nwcsequence{1} = nt_seq(2:end-1);
        elseif strcmp(loop_type,'IL')
            File.AllLoops_table.seq   = [nt_seq(1:breaks)       '*' nt_seq(breaks+1:end)];
            File.AllLoops_table.r_seq = [nt_seq(breaks+1:end)   '*' nt_seq(1:breaks)];
            File.AllLoops_table.nwc   = [nt_seq(2:breaks-1)     '*' nt_seq(breaks+2:end-1)];
            File.AllLoops_table.r_nwc = [nt_seq(breaks+2:end-1) '*' nt_seq(2:breaks-1)];
            File.AllLoops_table.sequence{1} = [nt_seq(1:breaks)       '*' nt_seq(breaks+1:end)];
            File.AllLoops_table.sequence{2} = [nt_seq(breaks+1:end)   '*' nt_seq(1:breaks)];
            File.AllLoops_table.nwcsequence{1} = [nt_seq(2:breaks-1)     '*' nt_seq(breaks+2:end-1)];
            File.AllLoops_table.nwcsequence{2} = [nt_seq(breaks+2:end-1) '*' nt_seq(2:breaks-1)];
        elseif strcmp(loop_type(1),'J')
            rangenumbers = [1:length(ranges) 1:length(ranges)];  % easier to cycle
            for rotation = 1:length(ranges),
                sequence{rotation} = nt_seq(ranges{rotation});
                nwcsequence{rotation} = nt_seq(nwcranges{rotation});
                for kk = 1:length(ranges),
                    sequence{rotation} = [sequence{rotation} '*' nt_seq(ranges{rangenumbers(rotation+kk)})];
                    nwcsequence{rotation} = [nwcsequence{rotation} '*' nt_seq(nwcranges{rangenumbers(rotation+kk)})];
                end
                File.AllLoops_table.sequence{rotation} = sequence{rotation};
                File.AllLoops_table.nwcsequence{rotation} = nwcsequence{rotation};
            end
            File.AllLoops_table.seq   = sequence{1};
            File.AllLoops_table.r_seq = sequence{2};       % for backward compatibility
            File.AllLoops_table.nwc   = nwcsequence{1};
            File.AllLoops_table.r_nwc = nwcsequence{2};
        end

        if i == 1
            Loops(1:length(Search.Candidates(:,1))) = File;
        else
            Loops(i) = File; %#ok<AGROW>
        end

    end
end

%==========================================================================
% Run HL or JL searches depending on the query passed.
%==========================================================================
function [Loops] = aRunHLorJLSearches(Query,pdb_id)

    Query.SearchFiles = pdb_id;
    Search = aFR3DSearch(Query);
    Loops = aProcessSearchData(Search);

end

%==========================================================================
% Search processing common to all types of searches.
%==========================================================================
function [Loops] = aProcessSearchData(Search)

    if isfield(Search,'Candidates') && ~isempty(Search.Candidates)
        Search = aFilterDuplicates(Search);
        Search.File = zAddNTData(Search.Query.SearchFiles); % load temporarily
        Loops = aSearchToLoops(Search);
    else
        Loops = 0;
    end

end

%==========================================================================
% Internal loop searches.
%==========================================================================
function [Loops] = aCombinedInternalLoopSearch(pdb_id)

    Query1 = aSetUpAlphaLoopQuery();
    Query2 = aSetUpStandardInternalLoopQuery();

    fprintf('%s\n',pdb_id);
    Query1.SearchFiles = pdb_id;
    S1 = aFR3DSearch(Query1);

    Query2.SearchFiles = pdb_id;
    S2 = aFR3DSearch(Query2);

    Search = aMergeSearches(S1,S2);

    Loops = aProcessSearchData(Search);

end

function [Query] = aSetUpAlphaLoopQuery()

    Query = aSetUpGenericQuery();
    Query.Name = 'AlphaLoopSearch';

    Query.NumNT = 4;

    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Edges{1,4} = 'cWW';
    Query.Edges{2,3} = 'cWW';
    Query.Edges{3,4} = 'flank';

    Query.Diff  = cell(Query.NumNT,Query.NumNT);
    Query.Diff{2,1}  = '=1';

    Query.Diagonal(1:Query.NumNT) = {'N'};

end

function [Query] = aSetUpStandardInternalLoopQuery()

    Query = aSetUpGenericQuery();
    Query.Name = 'StandardIL';

    Query.NumNT = 4;

    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Edges{1,4} = 'cWW';
    Query.Edges{1,2} = 'flankSS';
    Query.Edges{2,3} = 'cWW';
    Query.Edges{3,4} = 'flank';

    Query.Diff = cell(Query.NumNT,Query.NumNT);

    Query.Diagonal(1:Query.NumNT) = {'N'};

end

function [Search] = aMergeSearches(S1,S2)

    if ~isfield(S1,'Candidates')
        S1.Candidates = [];
        S1.File = [];
    end
    if ~isfield(S2,'Candidates')
        S2.Candidates = [];
        S2.File = [];
    end

    Search = S1;
    Search.Candidates = [S1.Candidates; S2.Candidates];
    if isfield(S1,'File')
        Search.File = S1.File;
    elseif isfield(S2,'File')
        Search.File = S2.File;
    end

end

%==========================================================================
% Generic hairpin search.
%==========================================================================
function [Loops] = aFlankHairpinSearch(pdb_id)

    Query = aSetUpGenericQuery();
    Query.Name = 'Hairpins';

    Query.NumNT = 2;
    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Edges{1,2} = 'flank cWW';

    Query.Diff = cell(Query.NumNT,Query.NumNT);
    Query.Diff{2,1} = '>';

    Query.Diagonal = {'N','N'};

    Loops = aRunHLorJLSearches(Query,pdb_id);

end

%==========================================================================
% Generic junction searches.
%==========================================================================
function [Loops] = aFlankJunctionSearch(pdb_id,loopType)

    Query = aSetUpGenericQuery();
    Query.Name = 'Junctions';

    Query.NumNT = 2 * str2double(loopType(2)); % J3, J4 etc
    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Diff  = cell(Query.NumNT,Query.NumNT);

    for i = 1:Query.NumNT-1
        if mod(i,2) == 0
            Query.Edges{i,i+1} = 'flank';     % may only find strands with non-empty interior
%             Query.Diff{i+1,i}  = '>';
        else
            Query.Edges{i,i+1} = 'cWW';
        end
    end
    Query.Edges{1,end} = 'flank';
    Query.Diagonal(1:Query.NumNT) = {'N'};

    Loops = aRunHLorJLSearches(Query,pdb_id);

%     aSaveAndCleanUp(filesToSearch,{Query.Name});

end

%==========================================================================
% Auxiliary functions.
%==========================================================================
function [Query] = aSetUpGenericQuery()

    Query.SaveDir = pwd;
    Query.Geometric = 0;
    Query.ExcludeOverlap = 0;
    Query.Number = 1;
    Query.Verbose = 0;

end

% aGetRange sorts the indices in the given row and interprets them as
% starting and ending ranges or strands of nucleotides, because we are
% using the result of a search for flanking cWW basepairs
% The result, range, will be all indices of nucleotides on all strands.
function [range,breaks,ranges,nwcranges] = aGetRange(row)

    % accepts Search.Candidates(i,:) as an input
    row = sort(row(1:end-1)); % discard the pdbfile id
    c = length(row);

    if c == 2 % hairpins
        range = row(1):row(2);
        breaks = [];
        ranges{1} = row(1):row(2);
        nwcranges{1} = (row(1)+1):(row(2)-1);
    elseif c == 4 % internal loops
        range = [row(1):row(2) row(3):row(4)];
        breaks = row(2);
        ranges{1} = row(1):row(2);
        ranges{2} = row(3):row(4);
        nwcranges{1} = (row(1)+1):(row(2)-1);
        nwcranges{2} = (row(3)+1):(row(4)-1);
    elseif mod(c,2) == 0 % junctions
        range = [];
        breaks = [];
        for q = 1:2:(c-1)
            range = [range row(q):row(q+1)];  %#ok<AGROW>
            breaks = [breaks row(q+1)];
            ranges{(q+1)/2} = row(q):row(q+1);
            nwcranges{(q+1)/2} = (row(q)+1):(row(q+1)-1);
        end
    else
        error('Wrong number of nucleotides');
    end

end

function [] = aCleanUp()

    filesToDelete = {'Hairpins','StandardIL','Junctions','AlphaLoopSearch'};
    for i = 1:length(filesToDelete)
        filename = [filesToDelete{i} '.mat'];
        if exist(filename,'file')
            delete(filename);
        end
    end

end

function [Search] = aFilterDuplicates(Search)

    ToSort = Search.Candidates;
    for i = 1:length(ToSort(:,1))
        ToSort(i,1:end-1) = sort(ToSort(i,1:end-1));
    end
    Search.Candidates = unique(ToSort,'rows');

end
