function [disc] = pairwiseSearch(file1, file2)

    % defer loading loop Precomputed data if not necessary
    if ~ischar(file1)        
        [File1, file1] = getNameAndData(file1);
    end
    if ~ischar(file2)
        [File2, file2] = getNameAndData(file2);        
    end
    
    result = getSearchAddress(file1, file2);

    % Parameters structure
    P.Discrepancy   = 1;
    P.Subdir        = fullfile(getSearchFolder, file1, '');    
    P.no_candidates = [P.Subdir filesep 'No_candidates.txt'];
    P.maxNtToSearch = 25;
    if ~exist(P.Subdir,'dir'), mkdir(P.Subdir); end
    
    % look up existing results
    if exist(result,'file')
        try 
            load(result);
            disc = Search.Discrepancy(1);
            return;
        catch
            fprintf('Corrupted file %s\n',result);
        end
    else
        P.no_candidates = fullfile(getSearchFolder, file1, 'No_candidates.txt');
        if exist(P.no_candidates, 'file')
        	fid = fopen(P.no_candidates, 'r');
        	no_candidates = textscan(fid, '%s');
        	fclose(fid);
            if ~isempty(find(ismember(no_candidates{1}, file2), 1))
                disc = Inf;
                return;
            end            
        end
    end
    
    % if file1 and file2 are never compared before loop ids
    if ~exist('File1', 'var')
        [File1, file1] = getNameAndData(file1);        
    end
    if ~exist('File2', 'var')
        [File2, file2] = getNameAndData(file2);        
    end
 
    % do not search some large loops because they cause Matlab to run out of memory
    % TODO: replace this temporary fix with a more general solution  
    if strcmp(file1, 'IL_1NBS_003') || strcmp(file1, 'IL_1NBS_010') || strcmp(file2, 'IL_1NBS_003') || strcmp(file2, 'IL_1NBS_010')        
	disc = Inf;
	addToNoCandidatesFile(file2, P);
	return;
    end

    % don't analyze huge spurious loops
    if File1.NumNT > P.maxNtToSearch || File2.NumNT > P.maxNtToSearch
        fprintf('Large loop: %s vs %s, %i vs %i\n', file1, file2, File1.NumNT, File2.NumNT);
        if file1 == file2
            P.Discrepancy = 0.1; % can search in itself with very low disc
        else
            addToNoCandidatesFile(file2, P);
            disc = Inf;
            return;
        end
    end
    
    % don't search large in small
    if ~isSizeCompatible(File1, File2)
        fprintf('Query %s is larger than target %s\n', file1, file2);
        addToNoCandidatesFile(file2, P);
        disc = Inf;
        return;
    end    
    
    % try searching
    Search = struct;
    
    if aSearchFlankingBases(File1,File2,P) == 1   
    
        [Query,S] = aConstructPairwiseSearch(File1, File2, P);

        % S.File(1) - query file, S.File(2) - target file.        
        if Query.NumNT <= S.File(2).NumNT
            Search = aFR3DSearch(Query,S);
        end
    end

    if ~isfield(Search,'Candidates') || isempty(Search.Candidates)        
        addToNoCandidatesFile(file2, P);        
        disc = Inf;
    else
        disc = min(Search.Discrepancy);
    end

end

function [searchPossible] = isSizeCompatible(File1, File2)

    if strcmp(File1.Filename, File2.Filename)
        searchPossible = 1;
        return;
    end

    L1 = File1.NumNT;
    B1 = length(aDetectBulgedBases(File1));
    queryEffectiveLength = L1 - B1;

    % check if the query without bulges is <= than the entire target loop
    if queryEffectiveLength <= File2.NumNT
        searchPossible = 1;
    else
        searchPossible = 0;
    end    

end

function [] = addToNoCandidatesFile(loop_id, P)

    fid = fopen(P.no_candidates, 'a');
    fprintf(fid,'%s\n', loop_id);
    fclose(fid); 

end

function [F, f] = getNameAndData(input_entity)

    if isstruct(input_entity) % assume File structures
        F = input_entity;
        f = input_entity.Filename;
    else % assume loop id
        load(getPrecomputedDataAddress(input_entity));
        F = File;
        f = File.Filename;
    end

end

function [matched] = aSearchFlankingBases(File1, File2, P)

    matched = 0;

    % skip hairpins
    if strcmp(File1.Filename(1:2),'HL')
        matched = 1;
        return;                
    end
    
    % don't search in the same file
    if File1.Filename == File2.Filename
        matched = 1;
        return;
    end        

    F = File1;
    [S.File(1), Indices] = leaveOnlyFlankingBases(File1);            
    S.File(2)            = leaveOnlyFlankingBases(File2);
    S.QIndex = [1 2];
    
    Query.Geometric      = 1;
    Query.ExcludeOverlap = 1;    
    Query.SaveDir        = pwd;
    Query.SearchFiles    = File2.Filename;
    Query.Filename       = F.Filename;                
    Query.ChainList      = {F.NT(Indices).Chain}; 
    Query.Name           = [File1.Filename '_' File2.Filename];        
    Query.NumNT          = length(Indices);
    Query.NTList         = {F.NT(Indices).Number};
    Query.NT             = F.NT(Indices);
    
    Query.Diagonal(1:Query.NumNT) = {'N'};
    Query.Edges = cell(Query.NumNT, Query.NumNT);
    Query.Edges{1,4} = 'cWW';
    Query.Edges{2,3} = 'cWW';
    
    Query.DiscCutoff = P.Discrepancy;                
    Query.RelCutoff  = Query.DiscCutoff;           
    
    Search = aFR3DSearch(Query,S);
    
    if isfield(Search,'Candidates') && ~isempty(Search.Candidates)
        matched = 1;
        if ismac
            unix(sprintf('rm %s', fullfile(pwd, [Query.Name '.mat'])));
        else
            delete([Query.Name '.mat']); % delete small search        
        end
    end            
    
end

function [File, Indices] = leaveOnlyFlankingBases(File)
    chbr = File.chain_breaks;
    Indices = [1 chbr chbr+1 length(File.NT)];    
    fn = fieldnames(File);
    for j = 1:length(fn)
        [r,c] = size(File.(fn{j}));
        if r == c && r == File.NumNT
            File.(fn{j}) = File.(fn{j})(Indices, Indices);
        end            
    end
    File.NT = File.NT(Indices);
    File.NumNT = 4;        
end

function [Query,S] = aConstructPairwiseSearch(File1, File2, P)
    
    F = File1;
    S.File = [File1 File2];
    S.QIndex = [1 2];

    if strcmp(F.Filename(1:2),'IL')
        bulges = aDetectBulgedBases(F);
        Indices = setdiff(1:length(F.NT),bulges);
    else
        Indices = 1:length(F.NT);
    end
                   
    Query.Geometric      = 1;
    Query.ExcludeOverlap = 1;    
    Query.SaveDir        = P.Subdir;
    Query.SearchFiles    = File2.Filename;
    Query.Filename       = F.Filename;    
    Query.ChainList      = {F.NT(Indices).Chain}; 
    Query.Name           = [File1.Filename '_' File2.Filename];        
    Query.Number         = 1;
    Query.NumNT          = length(Indices);
    Query.NTList         = {F.NT(Indices).Number};
    Query.NT             = F.NT(Indices);
    Query.IndicesManual  = Indices;
    
    Query.Diagonal(1:Query.NumNT) = {'N'};
    Query.Edges = cell(Query.NumNT,Query.NumNT);    
    
    Query.DiscCutoff = P.Discrepancy;                
    Query.RelCutoff = Query.DiscCutoff;           
    
    for i=1:Query.NumNT-1
        Query.Diff{i+1,i} = '>';       
    end                   

    if strcmp(File1.Filename(1:2), 'IL')
        Query.Edges{1,Query.NumNT} = 'cWW';        
        chainbreak = find(Indices==File1.chain_breaks);        
        Query.Diff{chainbreak+1,chainbreak} = '';    
        Query.Edges{chainbreak,chainbreak+1} = 'cWW';
        
        if File1.Flank(1,File1.chain_breaks) == 1
            Query.Edges{1, find(Indices==File1.chain_breaks)} = 'flankSS';
        end

        if File1.Flank(File1.chain_breaks+1, end) == 1
            Query.Edges{find(Indices==File1.chain_breaks+1), end} = 'flankSS';
        end
    else % HL
        Query.Edges{1, Query.NumNT} = 'cWW flankSS';
    end

end
