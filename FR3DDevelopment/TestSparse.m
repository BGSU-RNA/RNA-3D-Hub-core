
N = 5000;

s = round(rand(N,1)*100);

i = ceil(rand(N,1)*N);
j = ceil(rand(N,1)*N);

A = sparse(i,j,s,N,N);

i = ceil(rand(N,1)*N);
j = ceil(rand(N,1)*N);

B = sparse(i,j,s,N,N);

i = ceil(rand(N,1)*N);
j = ceil(rand(N,1)*N);

C = sparse(i,j,s,N,N);

tic
D = A + B;
toc
tic
E = A * B;
toc
tic
F = A .* B;
toc
tic
W = zSparseRange(C,0.7,0.9);
toc
tic
V = C .* (C <= 0.9) .* (C >= 0.7);
toc


H = zSparseValues(C,[1:30]);

