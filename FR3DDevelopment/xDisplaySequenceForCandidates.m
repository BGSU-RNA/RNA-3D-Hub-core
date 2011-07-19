% xDisplaySequenceForCandidates loads all available alignment information and summarizes the sequence variability for each candidate in a search

% We need a way to indicate which alignment is desired for each structure,
% since there may be several alignments.

function [File] = xDisplaySequenceForCandidates(File,Search,NofI)

MaxVariants = Inf;
IncludeGaps = 1;                             % include insertions after
                                             % nucleotides of interest?

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

Cand  = Search.Candidates;
Query = Search.Query;
[L,N] = size(Cand);
N = N - 1;                                       % number of nucleotides

if nargin < 3,
  NofI = 1:N;                                    % nucleotides of interest.
end

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

if ~isfield(Search,'Filenames'),
  for f = 1:length(Search.File),
    Search.Filenames{f} = Search.File(f).Filename;
  end
end

FN = unique(Search.Candidates(:,N+1));           % files with candidates

if ~isempty(File),
  [File,FIndex] = zAddNTData(Search.Filenames(FN),0,File,1);
else
  File = zAddNTData(Search.Filenames(FN),0,[],1);
end

% ------------------------------------------------ Attach alignment if needed

for f = 1:length(FIndex),
  if ~isfield(File(FIndex(f)).NT(1),'FASTA'),
    File(FIndex(f)) = zAttachAlignment(File(FIndex(f)),1);          
  end
end

% ------------------------------------------------ Loop through candidates

for c = 1:L,                                     % loop through candidates
 f = find(Search.Candidates(c,N+1)==FN);         % new file number
 f = Search.Candidates(c,N+1);                   % file number

 s = [];

 if isfield(File(f).NT(Cand(c,1)),'FASTA'),  % File(f) has sequence info

  for j = 1:N,                               % get FASTA column lengths
    s(j) = length(File(f).NT(Cand(c,j)).FASTA);
  end

  if min(s) == max(s) && min(s) > 0,         % lengths must be the same
    seq = [];
    for j = NofI,                            % paste together FASTA columns
      seq = [seq File(f).NT(Cand(c,j)).FASTA]; % append one column
      if IncludeGaps > 0 && j < NofI(end),
        seq = [seq File(f).NT(Cand(c,j)).GapsAfter]; % append more columns
      end
    end
    
    fprintf('Candidate %d is',c);
    for j = 1:N,
      fprintf(' %s%4s', File(f).NT(Cand(c,j)).Base,File(f).NT(Cand(c,j)).Number);
    end
    fprintf(' from %s chain %s\n', File(f).Filename, File(f).NT(Cand(c,1)).Chain);
    fprintf('Candidate %d has relevant sequence %s\n', c, cat(2,File(f).NT(Cand(c,NofI)).Base));

    fprintf('Summary of sequence variants found in alignment:\n');

    [b,t] = zUniqueRows(seq);

    L = min(MaxVariants,length(t));

    for r = 1:L,
      fprintf('Sequence %s occurs %4d times\n', b(r,:), t(r));
    end

    if L < length(t),
      fprintf('The remaining variants account for %8.2f%% of the sequences\n', 100*(1-sum(t(1:L))/sum(t)));
    end

    fprintf('\n');
  end
 end
end

