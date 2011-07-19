% zInferDesiredSSInSearch(Search) attempts to infer which "single-stranded" regions between nucleotides in the search should be kept, and stores these in the variable Between

function [Between] = zInferDesiredSSInSearch(Search)

[L,N] = size(Candidates);
N = N - 1;

Between = zeros(N,N);                            % default is no SS regions

% ---------------------------- Use MaxDiff constraints to infer

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(p,p),1);
else
  MaxDiff = ones(1,N-1);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
end

% ---------------------------- Check about Truncate

if ~isfield(Search,'Truncate'),
  Search.Truncate = [];
else
  maxinsert = 0*maxinsert;
end

% ---------------------------- Assign Between values

for n = 1:(N-1),                                % print alignment for cand
  if (maxinsert(n) < 4),
    Between(n,n+1) = 1;                         % pretty much contiguous
  end
  if (MaxDiff(n) < Inf),
    Between(n,n+1) = 1;                         % pretty much contiguous
  end
  if any(n+1 == Truncate),
    Between(n,n+1) = 0;
  end
end
