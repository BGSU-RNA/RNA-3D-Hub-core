% zDistance(A,B) finds the Euclidean distances between the rows of A and of B

function [D] = zDistance(A,B)

[M,t] = size(A);

if nargin == 2,
  [N,s] = size(B);
else
  N = 0;
  s = 0;
end

D = [];

if nargin < 2 && M < 150,

  G = A * A';                        % inner products of rows of A

  X = diag(G) * ones(1,M);           % M by M matrix

  D = sqrt(X + X' - 2*G);            % |u-v| = sqrt(|u|^2 + |v|^2 - 2 u . v)

elseif nargin < 2,

  G = A * A';                        % inner products of rows of A

  X = repmat(diag(G),1,M);           % repeat a in each column

  D = sqrt(X + X' - 2*G);            % |u-v| = sqrt(|u|^2 + |v|^2 - 2 u . v)

elseif M < 150 && s == t,

  a = sum((A.*A)');                   % sum of squares of each row of A
  b = sum((B.*B)');                   % sum of squares of each row of B

  X = a' * ones(1,N);                 % M by N matrix
  Y = ones(M,1) * b;                  % M by N matrix

  D = sqrt(X + Y - 2*A*B');           % |u-v| = sqrt(|u|^2 + |v|^2 - 2 u . v)

elseif s == t,

  [N,s] = size(B);

  a = sum((A.*A)');                   % sum of squares of each row of A
  b = sum((B.*B)');                   % sum of squares of each row of B

  X = repmat(b,M,1);                  % repeat b in each row
  Y = repmat(a',1,N);                 % repeat a in each column

  D = sqrt(X + Y - 2*A*B');           % |u-v| = |u|^2 + |v|^2 - 2 u . v

else

  D = [];
  fprintf('zDistance: Matrix sizes are not compatible\n');

end
