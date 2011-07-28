
K = 180;
W = 4000;

tic
for n = 1:W,
  A = rand(K,3);
  B = rand(K,3);
  D = zDistance(A,B);
end
toc


%break

tic
for n = 1:W,
  A = rand(K,3);
  B = rand(K,3);
  D = zDistanceOld(A,B);
end
toc

break

tic
for n = 1:W,
  A = rand(K,3);
  B = rand(K,3);
  t = cputime;
  D = zDistanceOld(A,B);
  elapse(n) = cputime - t;
end
toc


% ---------- Test that they give the same distance matrix

K = 18;
W = 40;

clear dif

for n = 1:W,
  A = rand(K,3);
  B = rand(K,3);
  D = zDistance(A,A);
  E = zDistanceOld(A,A);

  dif(n) = max(max(D-E));
end

sort(dif)
