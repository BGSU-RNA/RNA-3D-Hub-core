function [] = aEstimateBasepairVsTripleFrequency(runSearches)

if nargin < 1
    runSearches = 0;
end

% basepairs = {'tSS'};
basepairs = {'cWS', 'tWS', 'cSS', 'tSS', 'tHS'};

fid = fopen(['Frequency' filesep 'BasepairsVsTriples.html'],'w');
fprintf(fid,'<html><head><style type="text/css">');
fprintf(fid,'table {border-collapse:collapse;border: 1px solid;margin-left:auto;margin-right:auto;}');
fprintf(fid,'td {border: 1px solid;text-align:right;}');
fprintf(fid,'th {border: 1px solid;text-align:center;}');
fprintf(fid,'</style></head><body>');

for i = 1:length(basepairs)

    clear('Query','Search','SearchTp','SearchBp');
    if runSearches == 1
        Query = aCreateBasepairQuery(basepairs{i});
        SearchBp = aFR3DSearch(Query);
        clear('Query','Search');
        Query = aCreateTripleQuery(basepairs{i});
        SearchTp = aFR3DSearch(Query);
    else
        basepair = ['Frequency' filesep basepairs{i} filesep basepairs{i} '_basepair.mat'];
        load(basepair);
        SearchBp = Search;

        triple = ['Frequency' filesep basepairs{i} filesep basepairs{i} '_triple.mat'];
        load(triple);
        SearchTp = Search;
    end

    [freqBp,l1] = aBreakDownByNucleotide(SearchBp);    
    SearchTp = aRemoveMultipleTriples(SearchTp);
    [freqTp,l2] = aBreakDownByNucleotide(SearchTp);

    freq = freqTp ./ freqBp; 

    disp(freq);
    fprintf('\n\n');
%     dlmwrite(['Frequency' filesep basepairs{i} filesep 'freq.txt'], freq, 'delimiter', '\t', 'precision', 4);
%     dlmwrite(['Frequency' filesep basepairs{i} filesep 'freqTp.txt'], freqTp, 'delimiter', '\t', 'precision', 4);
%     dlmwrite(['Frequency' filesep basepairs{i} filesep 'freqBp.txt'], freqBp, 'delimiter', '\t', 'precision', 4);

    aWriteFrequenciesToFile(fid,freq,freqBp,basepairs{i});
    
end

fprintf(fid,'</body></html>');
fclose(fid);
keyboard;

end

function [Search] = aRemoveMultipleTriples(Search)

    Search.OldCandidates = Search.Candidates;
    pdbcolumn = length(Search.Candidates(1,:));
    Search.Candidates = sortrows(Search.Candidates,[pdbcolumn 1]);
    i = 1;
    removed = 0;
    while i < length(Search.Candidates(:,1))
        if Search.Candidates(i,[1 2 4]) == Search.Candidates(i+1,[1 2 4])
            Search.Candidates(i+1,:) = [];
            removed = removed + 1;
        else
            i = i + 1;
        end
    end
    fprintf('Removed %i entries\n',removed);
    
end

function [freq,l] = aBreakDownByNucleotide(Search)

    freq = zeros(4,4);
    l = [];
    for index = 1:length(Search.Candidates(:,1))
        pdbfile = Search.Candidates(index,end);        
        i = Search.File(pdbfile).NT(Search.Candidates(index,1)).Code;
        j = Search.File(pdbfile).NT(Search.Candidates(index,2)).Code;
        if i == 1 && j == 2
            l = [l index];
        end
        freq(i,j) = freq(i,j) + 1;
    end

end

function [i,j] = aFirstTwoNucleotides(Search, index) %#ok<DEFNU>

    pdbfile = Search.Candidates(index,end);
    nt(1) = Search.File(pdbfile).NT(Search.Candidates(index,1)).Base;
    nt(2) = Search.File(pdbfile).NT(Search.Candidates(index,2)).Base;

    for q = 1:length(nt)
        switch nt(q)
            case 'A' 
                number(q) = 1;
            case 'C' 
                number(q) = 2;
            case 'G' 
                number(q) = 3;
            case 'U' 
                number(q) = 4;           
        end    
    end

    i = number(1);
    j = number(2);

end

function [Query] = aCreateBasepairQuery(basepair)

    Location = ['Frequency' filesep basepair];
    if ~exist(Location, 'dir')
        mkdir(Location);
    end     

    Query.SaveDir = Location;
    Query.Geometric = 0;
    Query.ExcludeOverlap = 1;
    Query.Name = [basepair '_basepair'];
    Query.NumNT = 2;
    Query.Number = 1;
    Query.Edges = cell(2,2);
    Query.Edges{1,2} = basepair;
    Query.Diff = cell(2,1);
    Query.Diagonal = {'N','N'};
    Query.SearchFiles{1} = 'NAR_NR_list.pdb';    
%     Query.SearchFiles{1} = 'NR_list_2009-05-14.pdb';
    % Query.SearchFiles{1} = '1J5E';
end

function [Query] = aCreateTripleQuery(basepair)

    Location = ['Frequency' filesep basepair];
    if ~exist(Location, 'dir')
        mkdir(Location);
    end     

    Query.SaveDir = Location;
    Query.Geometric = 0;
    Query.ExcludeOverlap = 1;
    Query.Name = [basepair '_triple'];
    Query.NumNT = 3;
    Query.Number = 1;
    Query.Edges = cell(3,3);
    Query.Edges{1,2} = basepair;
    Query.Edges{2,3} = 'pair';
    Query.Diff = cell(3,2);
    Query.Diagonal = {'N','N','N'};
    Query.SearchFiles{1} = 'NAR_NR_list.pdb';        
%     Query.SearchFiles{1} = 'NR_list_2009-05-14.pdb';
    % Query.SearchFiles{1} = '1J5E';
end

function [] = aWriteFrequenciesToFile(fid,freq,freqBp,basepair)

    letters = 'ACGU';
    fprintf(fid,'<table><caption>%s</caption><tr><td>&nbsp</td><th colspan="2">A</th>',basepair);
    fprintf(fid,'<th colspan="2">C</th><th colspan="2">G</th><th colspan="2">U</th></tr>');
    fprintf(fid,'<tr><td>&nbsp</td><td>Triple/Basepair</td><td>Basepair</td><td>Triple/Basepair</td><td>Basepair</td>');
    fprintf(fid,'<td>Triple/Basepair</td><td>Basepair</td><td>Triple/Basepair</td><td>Basepair</td></tr>');
    for j = 1:4
        for k = 1:5
            if k == 1
                fprintf(fid,'<tr><th>%s</th>',letters(j));
            else
                if isnan(freq(j,k-1))
                    fprintf(fid,'<td>0.0</td>');                    
                else
                    fprintf(fid,'<td>%.1f</td>',freq(j,k-1)*100);
                end
                fprintf(fid,'<td>%i</td>',freqBp(j,k-1));                
            end
        end
        fprintf(fid,'</tr>');
    end
    fprintf(fid,'</table><br><br>');

end

