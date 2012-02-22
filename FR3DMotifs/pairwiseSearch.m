function [disc] = pairwiseSearch(file1, file2)

    [File1, file1] = getNameAndData(file1);
    [File2, file2] = getNameAndData(file2);    

    name = getSearchAddress(file1, file2);
    
    if exist(name,'file')
        try 
            load(name);
            disc = Search.Discrepancy(1);
            return;
        catch
            fprintf('Corrupted file %s\n',name);
        end
    end

    % Parameters structure
    P.Discrepancy  = 1;
    P.Subdir       = fullfile(getSearchFolder, file1, '');    
    Search = struct;
    
    if aSearchFlankingBases(File1,File2,P) == 1   
    
        [Query,S] = aConstructPairwiseSearch(File1, File2, P);

        % S.File(1) - query file, S.File(2) - target file.        
        if Query.NumNT <= S.File(2).NumNT
            Search = aFR3DSearch(Query,S);
        end
    end

    if ~isfield(Search,'Candidates') || isempty(Search.Candidates)

        if ~exist(P.Subdir,'dir'), mkdir(P.Subdir); end
        
        fid = fopen(fullfile(P.Subdir, 'No_candidates.txt'),'a');        
        fprintf(fid,'%s\n',file2);
        fclose(fid); 
        disc = Inf;
    else
        disc = min(Search.Discrepancy);
    end

end

function [F, f] = getNameAndData(input_entity)

    if isstruct(input_entity)
        F = input_entity;
        f = input_entity.Filename;
    else
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
        delete([Query.Name '.mat']); % delete small search        
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

    Query.Edges{1,Query.NumNT} = 'cWW';
    
    chainbreak = find(Indices==File1.chain_breaks);
    if ~isempty(chainbreak) %IL
        Query.Diff{chainbreak+1,chainbreak} = '';    
        Query.Edges{chainbreak,chainbreak+1} = 'cWW';
        
        if File1.Flank(1,File1.chain_breaks) == 1
            Query.Edges{1, find(Indices==File1.chain_breaks)} = 'flankSS';
        end

        if File1.Flank(File1.chain_breaks+1, end) == 1
            Query.Edges{find(Indices==File1.chain_breaks+1), end} = 'flankSS';
        end
    end

end
