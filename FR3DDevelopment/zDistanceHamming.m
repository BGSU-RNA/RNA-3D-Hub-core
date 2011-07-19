% zDistance(A,B) finds the Euclidean distances between the rows of A and of B
% A and B need to have the same number of columns
% For faster performance, let A be the larger matrix, B be the smaller one

function [D] = zDistance(A,B)

[M,t] = size(A);

if nargin < 2,

  D = zDistanceHamming(A,A);

else

  [a,b] = size(A);
  [c,d] = size(B);

  if b ~= d,
    fprintf('zDistanceHamming: Rows are not the same lengths.\n');
  end

  D = zeros(a,c,'uint8');                   % D(i,j) is distance between
                                            % row i of A and row j of B

  for j = 1:c,                              % run through rows of B
    for k = 1:b,                            % run through columns
      D(:,j) = D(:,j) + uint8(A(:,k) ~= B(j,k)); % add 
    end
  end

end
