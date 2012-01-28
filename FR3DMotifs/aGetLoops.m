%==========================================================================
% Extracts all RNA 3D motifs of a specified loopType from a set of PDB files.
% Input: pdb id, motif loopType (IL, HL, J3, J4 etc).
% Returns loops as a structure array or 0 if no loops were found.
% All motifs, including noncontiguous and ncWW are stored. Postprocessing 
% should be done separately.
%==========================================================================
function [Loops, L, err_msg] = aGetLoops(pdb_id,loopType)

    global loop_type;

    try
        loop_type = upper(loopType);
        Loops     = 0;
        L         = 0;
        err_msg   = '';        
        
        if strcmp(loopType,'IL')
            Loops = aCombinedInternalLoopSearch(pdb_id);        
        elseif strcmp(loopType,'HL')
            Loops = aFlankHairpinSearch(pdb_id);
        elseif strfind(loopType,'J')
            Loops = aFlankJunctionSearch(pdb_id,loopType);
        end

        L = length(Loops);
        aCleanUp;

    catch err
        err_msg = sprintf('Error "%s" in aGetLoops on line %i (%s)\n', ...
                           err.message, err.stack.line, pdb_id);
        disp(err_msg);
    end
        
end

%==========================================================================
% Write out pdbs and small mat files for each loop.
%==========================================================================
function [Loops] = aWriteMatAndPdbFiles(Search)
        
    global loop_type;

    for i = 1:length(Search.Candidates(:,1))
        
        pdbfile = Search.Candidates(i,end);
        File    = Search.File(pdbfile);
        range   = aGetRange(Search.Candidates(i,:));
        
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
            ids{j} = [aGetNTId(File,j) ','];        
        end
        ids = sort(ids);
        File.AllLoops_table.full_id = [ids{:}];
        File.AllLoops_table.full_id = File.AllLoops_table.full_id(1:end-1);       
        
        if strcmp(loop_type,'IL')
            File.chain_breaks = find(range==Search.Candidates(i,2));
        elseif strcmp(loop_type,'J3')
            File.chain_breaks = find(range==Search.Candidates(i,2) | ...
                                     range==Search.Candidates(i,4));
        elseif strcmp(loop_type,'HL') % a hairpin
            File.chain_breaks = [];
        end        
        
        File.AllLoops_table.loop_name = aGetLoopName(File);        
        
        nt_seq   = [File.NT.Base];
        if strcmp(loop_type,'IL')
            breaks = File.chain_breaks;            
            File.AllLoops_table.seq   = [nt_seq(1:breaks)       '*' nt_seq(breaks+1:end)];
            File.AllLoops_table.r_seq = [nt_seq(breaks+1:end)   '*' nt_seq(1:breaks)];                   
            File.AllLoops_table.nwc   = [nt_seq(2:breaks-1)     '*' nt_seq(breaks+2:end-1)];
            File.AllLoops_table.r_nwc = [nt_seq(breaks+2:end-1) '*' nt_seq(2:breaks-1)];
        elseif strcmp(loop_type(1),'J')
            File.AllLoops_table.seq   = nt_seq;
            File.AllLoops_table.r_seq = '';
            File.AllLoops_table.nwc   = ''; %seq(2:end-1);
            File.AllLoops_table.r_nwc = '';            
        elseif strcmp(loop_type,'HL')
            File.AllLoops_table.seq   = nt_seq;
            File.AllLoops_table.r_seq = '';
            File.AllLoops_table.nwc   = File.AllLoops_table.seq(2:end-1);
            File.AllLoops_table.r_nwc = '';
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
        Loops = aWriteMatAndPdbFiles(Search);    
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
            Query.Edges{i,i+1} = 'flank';
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

function [range] = aGetRange(row)

    % accepts Search.Candidates(i,:) as an input
    c = length(row);
    row = sort(row(1:end-1)); % discard the pdbfile id
    
    if c == 3 % hairpins
        range = row(1):row(2);
    elseif c == 5 % internal loops
        range = [row(1):row(2) row(3):row(4)];
    elseif mod(c,2) == 1 % junctions
        range = [];
        for q = 1:2:c-1
            range = [range row(q):row(q+1)];  %#ok<AGROW>
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
