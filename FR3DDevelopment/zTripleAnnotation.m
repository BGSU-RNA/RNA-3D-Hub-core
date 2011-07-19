
if 0 > 1,
%load 2009-07-09_16_37_30-Base_triples.mat
load 2009-07-09_17_08_51-Base_triples.mat

File = zAddNTData(Search.Filenames);

xAnnotate
end

[L,N] = size(Search.Candidates);

N = N - 1;

r = 2;
clear T
clear Y

for c = 1:L,
  f = Search.Candidates(c,N+1);            % file number
  i = Search.Candidates(c,1);
  j = Search.Candidates(c,2);
  k = Search.Candidates(c,3);

  if fix(abs(File(f).Edge(i,j))) <= fix(abs(File(f).Edge(j,k))),
    T{r, 1} = r;
    T{r, 2} = File(f).NT(i).Base;
    T{r, 3} = File(f).NT(i).Number;
    T{r, 4 } = File(f).NT(j).Base;
    T{r, 5} = File(f).NT(j).Number;
    T{r, 6} = File(f).NT(k).Base;
    T{r, 7} = File(f).NT(k).Number;
    T{r, 8} = zEdgeText(File(f).Edge(i,j));
    T{r, 9} = zEdgeText(File(f).Edge(j,k));
    T{r,10} = zEdgeText(File(f).Edge(i,k));

    Y(r,1) = abs(File(f).Edge(i,j));
    Y(r,2) = abs(File(f).Edge(j,k));

    T{r,11} = File(f).Crossing(i,j);
    T{r,12} = File(f).Crossing(j,k);
    T{r,13} = File(f).Crossing(i,k);
    if ~isempty(File(f).Nucl(i).Motif),
      T{r,14} = strrep(File(f).Nucl(i).Motif(1).Name,'-','');
      Y(r,3) = File(f).Nucl(i).Motif(1).Index;
    end
    if ~isempty(File(f).Nucl(j).Motif),
      T{r,15} = strrep(File(f).Nucl(j).Motif(1).Name,'-','');
      Y(r,4) = File(f).Nucl(j).Motif(1).Index;
    end
    if ~isempty(File(f).Nucl(k).Motif),
      T{r,16} = strrep(File(f).Nucl(k).Motif(1).Name,'-','');
      Y(r,5) = File(f).Nucl(k).Motif(1).Index;
    end
    T{r,17} = File(f).Filename;
    r = r + 1;
  end
end

T{1,1} = 'Index';
T{1,2} = 'Base 1';
T{1,3} = 'Num 1';
T{1,4} = 'Base 2';
T{1,5} = 'Num 2';
T{1,6} = 'Base 3';
T{1,7} = 'Num 3';
T{1,8} = 'Pair 1-2';
T{1,9} = 'Pair 2-3';
T{1,10} = 'Pair 1-3';
T{1,11} = 'Range 1-2';
T{1,12} = 'Range 2-3';
T{1,13} = 'Range 1-3';
T{1,14} = 'SecStruct 1';
T{1,15} = 'SecStruct 2';
T{1,16} = 'SecStruct 3';
T{1,17} = 'PDB ID';

[a,b] = sortrows(Y,[1 2 3 4 5]);

T = T(b,:);

T = T(:,[1 17 2:16]);

xlswrite('TripleAnnotation.xls',T);
