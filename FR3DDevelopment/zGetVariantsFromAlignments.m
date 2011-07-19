% zGetVariantsFromAlignments(File,Search,Verbose,UseAlignments) 

function [File,seqs,codes,counts] = zGetVariantsFromAlignments(File,Search,Verbose,Truncate,UseAlignments)

if nargin < 3,
  Verbose = 0;
end

if nargin < 4,
  for f = 1:length(File),
    UseAlignments{f} = File(f).Filename;
  end
end

Cand = Search.Candidates;

[L,N] = size(Cand);                      % number of candidates
N = N - 1;                               % number of nucleotides

ff = unique(Cand(:,end));                % file numbers needed for candidates

Filenames(ff) = Search.Filenames(ff);    % file names needed

if ~isfield(File(ff(1)).NT(1),'FASTA')

  for f = 1:length(File),

File(f).Filename

    if any(ismember(UseAlignments,File(f).Filename)),

disp('Attaching');

      File(f) = zAttachAlignment(File(f),1);           % attach alignment data
    end
  end
end

AllLett = [];
clear Lett

for c = 1:L,                                % loop through candidates
  f = Cand(c,N+1);                          % file number
  Let = [];
  if isfield(File(f).NT(Cand(c,1)),'FASTA'),
    for n = 1:N,
      Let = [Let File(f).NT(Cand(c,n)).FASTA];
      Gaps{c,n} = File(f).NT(Cand(c,n)).GapsAfter; % gaps after this column
    end
  end

  AllLett = [AllLett; Let];                 % accumulate for later analysis

  % ----------------------------------------- analyze this instance

  [seqs,counts] = zUnique(Let);             % get unique sequences and counts

  codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');
  i = find(min(codes,[],2) > 0);            % identify sequences w/o non ACGU
  Let = seqs(i,:);

  Counts{c} = counts(i);
  Lett{c} = Let;

  if Verbose > 0,
    fprintf('File %s instance ', File(f).Filename);
    for n = 1:N,
      fprintf('%s%4s ', File(f).NT(Cand(c,n)).Base, File(f).NT(Cand(c,n)).Number);
    end
    fprintf('\n');

    for i = 1:length(Counts{c}),
      fprintf('  %s occurs %4d times\n', Lett{c}(i,:), Counts{c}(i));
    end

  end
end

% ----------------------------------- find unique sequences over all instances

[seqs,counts] = zUnique(AllLett);
codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');

% ----------------------------------- add in gaps except where Truncated

for n = 1:N,
  maxgap = 0;
  for c = 1:L,
    [Z,gapsize] = size(Gaps{c,n});
    maxgap = max(gapsize,maxgap);
  end
  for c = 1:L,
    [Z,gapsize] = size(Gaps{c,n});
    while gapsize < maxgap,
      [Z,gapsize] = size(Gaps{c,n});
    end
  end
end

% ------------------------------------------- remove non-ACGU letters

i = find(min(codes,[],2) > 0);
seqs = seqs(i,:);
codes = codes(i,:);
counts = counts(i);
[S,T] = size(seqs);                       % S is number of unique sequences

if Verbose > 0,
  fprintf('\n');
  fprintf('Summary of sequence variants\n');
  for i = 1:S,
    fprintf('  %s occurs %4d times\n', seqs(i,:), counts(i));
  end
end
