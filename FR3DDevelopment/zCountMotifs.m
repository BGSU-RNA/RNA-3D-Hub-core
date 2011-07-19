
% load 2008-06-06_13_44_49-Internal_loops_flanked_by_nested_cWW_pairs.mat
Filenames = Search.Filenames;
Verbose = 1;

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,2,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,2,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);

Cand = Search.Candidates;

[L,N] = size(Cand);
N = N - 1;                 % number of nucleotides

clear Text

r = 1;
while r <= L/2,
  m = min(Cand(2*r,[1 2]));          % the smaller 
  M = max(Cand(2*r,[1 2]));
  i = m:M;

  mm = min(Cand(2*r,[3 4]));
  MM = max(Cand(2*r,[3 4]));
  j = mm:MM;

  f = Cand(2*r,N+1);

  [a,b,c] = find(File(f).Edge(i,j));
  c = abs(c);
  c = c .* (c > 0) .* (c < 13);             % "true" basepairs
  c = c(find(c));
  c = sort(c);

  if length(c) == 2,
    [a,b,c] = find(File(f).Edge(i,j));
    c = abs(c);
    c = c .* (c > 100) .* (c < 113);       % "near" basepairs
    c = c(find(c));
    c = sort(c);
    c = [1; 1; c];
  end

  if length(c) > 0,
    Text{r} = '';
    for g = 1:length(c),
      Text{r} = [Text{r} sprintf(' %s', zEdgeText(c(g)))];
    end
  else
    Text{r} = ' No interaction';
  end

%  fprintf('%s\n',Text{r});


  if Verbose > 1,
    zShowInteractionTable(File(f),[mm:MM m:M]);

    VP.Sugar = 1;
    clf
    zDisplayNT(File(f), [mm:MM m:M],VP)
    disp(Text{r})
    drawnow
    pause
  end

  r = r + 1;

end

[t,p] = sort(Text);
Text = Text(p);

[b,i,j] = unique(Text);

% i points to the unique rows of Text
% j maps 1:L to 1:length(b)

fprintf('Of the %d internal loops identified, there were %d unique sequence signatures.\n', L, length(i));

for r = 1:length(i),
  fprintf('%5d instances of %s\n', length(find(j==r)), b{r});
end
