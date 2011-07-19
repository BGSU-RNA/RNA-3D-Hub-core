% xCountCandidateInteractions(Search) tallies up how many candidates have 0, 1, 2, ... basepairs, stacks, etc.

function [void] = xCountCandidateInteractions(Search)

[L,N] = size(Search.Candidates);

N = N - 1;

bpcount = zeros(1,N^2);

for c = 1:L,
  f = Search.Candidates(c,N+1);
  E = abs(fix(Search.File(f).Edge));

  e = E(Search.Candidates(c,1:N),Search.Candidates(c,1:N));  % interactions
  e = triu(e);                                % just take one set

  bp = sum(sum((e < 13) .* (e > 0)));

  bpcount(bp+1) = bpcount(bp+1) + 1;
end

bp = 1;
while bp <= length(bpcount) && sum(bpcount(bp:end)) > 0,
  fprintf('%4d candidates have %2d basepairs\n', bpcount(bp), bp-1);
  bp = bp + 1;
end

fprintf('Total of %d = %d candidates\n', L, sum(bpcount));
