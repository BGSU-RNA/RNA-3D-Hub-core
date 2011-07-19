% zUniqueRows determines the unique rows of M and returns counts of how many instances of each row are present.  The unique rows are returned in the matrix b, and the counts are returned in t. 

function [b,t,i] = zUniqueRows(M,z,Verbose)

if nargin < 2,
  z = 1:length(M(1,:));
end

if nargin < 3,
  Verbose = 0;
end

[N,j] = sortrows(M);                  % sorted rows with index j

t(1) = 1;
n    = N(1,z);
c    = 1;
b(c,:) = N(c,:);
i(c)   = j(1);

for r = 2:length(N(:,1)),
  if all(n == N(r,z)),                  % same as previous row
    t(c) = t(c) + 1;                    % add a count for it
  else
    n = N(r,z);                         % new current row
    c = c + 1;
    b(c,:) = N(r,:);
    t(c) = 1;                           % current count
    i(c) = j(r);
  end
end

[y,k] = sort(-t);

b = b(k,:);                             % most populous first
t = t(k);

% [b,i,j] = unique(M,'rows');

if Verbose > 0,
  if strcmp(class(M),'char'),
    for a = 1:length(t),
      fprintf('%s appears %4d times\n', b(i,:), t(i));
    end
  end
end

