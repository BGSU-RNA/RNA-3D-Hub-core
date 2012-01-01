%==========================================================================
% Extracts all RNA 3D motifs of a specified loopType from a set of PDB files.
% Input: folder, motif loopType (IL, HL, 3WJ, 4WJ etc), PDB files (optional).
% By default 1S72 is searched. PDB files can be specified as a cell array
% or as a list (must be located in /PDBFiles).
% All motifs, including noncontiguous and ncWW are stored. Postprocessing 
% should be done separately.
% Use J3, J4 etc to search for junctions
%==========================================================================
function [Loops, L, err_msg] = aGetLoops(pdb_id,loopType)

    global loop_type;
    try
        loop_type = loopType;        
        Loops = [];
        L = 0;
        err_msg = '';        
        
        if strcmp(loopType,'IL')
            Loops = aCombinedInternalLoopSearch(pdb_id);        
        elseif strcmp(loopType,'HL')
            Loops = aFlankHairpinSearch(pdb_id);
        elseif strfind(loopType,'J')
            Loops = aFlankJunctionSearch(pdb_id,loopType);
        end

        L = length(Loops);

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
            File.chain_breaks = figure_out_later;
        else % a hairpin
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
        else % HL
            File.AllLoops_table.seq   = nt_seq;
            File.AllLoops_table.r_seq = '';
            File.AllLoops_table.nwc   = File.AllLoops_table.seq(2:end-1);
            File.AllLoops_table.r_nwc = '';
        end
                       
        if i == 1
            Loops = File;
        else
            Loops = [Loops File];
        end
                
    end        
end

%==========================================================================
% Run HL or JL searches depending on the query passed.
%==========================================================================
function [] = aRunHLorJLSearches(Query,filesToSearch)

    for i = 1:length(filesToSearch)
        Query.SearchFiles = filesToSearch{i};
        Search = aFR3DSearch(Query);
        if isfield(Search,'Candidates') && ~isempty(Search.Candidates)
            aProcessSearchData(Search);
        end        
        aSaveAndCleanUp(filesToSearch{i},{Query.Name});
    end    
    
end

%==========================================================================
% Search processing common to all types of searches.
%==========================================================================
function [Loops] = aProcessSearchData(Search)

    Search = aFilterDuplicates(Search);             
    Search.File = zAddNTData(Search.Query.SearchFiles); % load temporarily      
    Loops = aWriteMatAndPdbFiles(Search);    

end

%==========================================================================
% Internal loop searches.
%==========================================================================
function [Loops] = aCombinedInternalLoopSearch(pdb_id)
    
    Loops = struct();

    Query1 = aSetUpAlphaLoopQuery();        
    Query2 = aSetUpSixNTLoopQuery();                
            
    fprintf('%s\n',pdb_id);
    Query1.SearchFiles = pdb_id;
    S1 = aFR3DSearch(Query1);        

    Query2.SearchFiles = pdb_id;
    S2 = aFR3DSearch(Query2);    

    if isfield(S2,'Candidates') && ~isempty(S2.Candidates)
        S2.Candidates = S2.Candidates(:,[1 3 4 6 7]);
    end    

    Search = aMergeSearches(S1,S2); 
    if isfield(Search,'Candidates') && ~isempty(Search.Candidates)
        Loops = aProcessSearchData(Search);
    end
            
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

function [Query] = aSetUpSixNTLoopQuery()

    Query = aSetUpGenericQuery();
    Query.Name = 'SixNtIL';

    Query.NumNT = 6;

    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Edges{1,3} = 'flank';	
    Query.Edges{1,6} = 'cWW';	
    Query.Edges{3,4} = 'cWW';		
    Query.Edges{4,6} = 'flank';		

    Query.Diff = cell(Query.NumNT,Query.NumNT);
    Query.Diff{2,1} = '=1 >';
    Query.Diff{3,2} = '>';
    Query.Diff{5,4} = '>';		
    Query.Diff{6,5} = '=1 >';	

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
function [] = aFlankHairpinSearch(filesToSearch)
    
    Query = aSetUpGenericQuery();    
    Query.Name = 'Hairpins';
    
    Query.NumNT = 2;
    Query.Edges = cell(Query.NumNT,Query.NumNT);
    Query.Edges{1,2} = 'flank cWW';

    Query.Diff = cell(Query.NumNT,Query.NumNT);
    Query.Diff{2,1} = '>';
    
    Query.Diagonal = {'N','N'};

    aRunHLorJLSearches(Query,filesToSearch); 
            
end

%==========================================================================
% Generic junction searches.
%==========================================================================
function [] = aFlankJunctionSearch(filesToSearch,loopType)

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
    
    aRunHLorJLSearches(Query,filesToSearch);  
    
    aSaveAndCleanUp(filesToSearch,{Query.Name});    
        
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

function [] = aSaveAndCleanUp(filesToDelete)
        
    for i = 1:length(filesToDelete)
        filename = [Param.location filesep filesToDelete{i} '.mat'];
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
