% xDistinctCandidates(Search) sorts the indices of each candidate in Search, identifies one instance of each, and returns a Search with just these candidates

function [Search,n] = xDistinctCandidates(Search,Verbose)

if nargin < 2,
  Verbose = 0;
end

if strcmp(class(Search),'char'),
  load(Search)
end

[L,N] = size(Search.Candidates);

N = N - 1;

C = sort(Search.Candidates(:,1:N),2);

C = [C Search.Candidates(:,N+1)];

[D,i] = unique(C,'rows');                     % identify the unique rows of C

n = length(i);

if Verbose > 0,
  fprintf('This search has %d distinct candidates\n', length(i));
end

Search.Candidates = Search.Candidates(i,:);  % keep one of each
if isfield(Search,'Discrepancy'),
  Search.Discrepancy = Search.Discrepancy(i);
end


