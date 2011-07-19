
c1 = 1;          % NTs to tabulate
c2 = 6;

c1 = 2;          % NTs to tabulate
c2 = 5;

Cand = Search.Candidates;

[L,N] = size(Cand);

N = N - 1;                     % number of query nucleotides

Count = zeros(4,4);

for c = 1:L,
  f = Cand(c,N+1);

  code1 = Search.File(f).NT(Cand(c,c1)).Code;
  code2 = Search.File(f).NT(Cand(c,c2)).Code;

  Count(code1,code2) = Count(code1,code2) + 1;  
end

Count