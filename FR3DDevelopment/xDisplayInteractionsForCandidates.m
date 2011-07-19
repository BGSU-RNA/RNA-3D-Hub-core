% xDisplaySequenceForCandidates loads all available alignment information and summarizes the sequence variability for each candidate in a search

function [File] = xDisplaySequenceForCandidates(File,Search,NofI)

if strcmp(class(Search),'double'),

  % ---------------------------------------- Load a saved search by number
  g = Search;

  switch g
  case 1,
    load LIB00005_IL_tSH-tHH-cSH-tWH-tHS_sarcin-ricin.mat
    NofI = 1:13;                      % nucleotides of interest
  case 2,
    load 2009-01-28_16_32_47-GC_1BPh_swap.mat
    NofI = [1 2];                     % nucleotides of interest for comparison
  case 3,
    load 2009-07-21_14_13_55-2avy_triples
    NofI = 1:3;
  case 4,
    load('LIB00005_IL_tSH-tHH-cSH-tWH-tHS_sarcin%3Aricin')
    NofI = 1:13;
  case 5,
    load('2009-08-03_17_01_39-Sarcin-Ricin_13_composite_non-redundant')
    NofI = 1:13;
  end

elseif strcmp(class(Search),'char'),

  load(Search);

end

% ------------------------------------------------ Preliminary calculations

Cand = Search.Candidates;
Query = Search.Query;
[L,N] = size(Cand);
N = N - 1;                                       % number of nucleotides

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
end

% ------------------------------------------------ Load file information

FN = unique(Search.Candidates(:,N+1));           % files with candidates

if ~isempty(File),
  [File,FIndex] = zAddNTData(Search.Filenames(FN),0,File,1);
  File = File(FIndex);
else
  File = zAddNTData(Search.Filenames(FN),0,[],1);
end

% ------------------------------------------------ Attach alignment if needed

for f = 1:length(File),
  if ~isfield(File(f).NT(1),'FASTA'),
%    File(f) = zAttachAlignment(File(f),1);          
  end
end

% ------------------------------------------------ Loop through candidates

for n = 1:N,                                     % loop through nucleotides
 f  = find(Search.Candidates(1,N+1)==FN);
 NT = File(f).NT(Search.Candidates(1,n));
 fprintf('Nucleotide %d of the query is %s%4s.\n', n, NT.Base, NT.Number);

 for c = 1:L,                                     % loop through candidates
  f = find(Search.Candidates(c,N+1)==FN);         % new file number

  s = [];

  i = Search.Candidates(c,1:N);                   % indices

  e = File(f).Edge(i(n),:);
  e(i) = zeros(1,N);                              % remove inter. with cand.

  j = find(e);                                    % interaction partners    

  NT = File(f).NT(i(n));

  fprintf('  Candidate %d from %s; this is %s%4s.\n', c, File(f).Filename, NT.Base, NT.Number);

  for k = 1:length(j),
    NT2 = File(f).NT(j(k));
    fprintf('    %4s with %s%4s\n', zEdgeText(e(j(k))), NT2.Base, NT2.Number);
  end
 end
end

