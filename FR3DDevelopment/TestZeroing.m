
if ~exist('File'),
  File = zAddNTData('1S72');
  c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
  DD = zMutualDistance(c,16); % compute distances < 16 Angstroms
end


E = File.Edge;

G = zSparseRange(abs(E),21,24);
D = DD;
nnz(D)

tic

[i,j] = find(G);
for k = 1:length(i),
  D(i(k),j(k)) = 0;
end
nnz(D)

toc

D = DD;

tic

D = D .* (G == 0);
nnz(D)

toc

D = DD;

tic

[i,j,d] = find(D);
[p,q] = find(G);

[c,k] = setdiff([i j], [p q], 'rows');

D = sparse(i(k),j(k),d(k));
nnz(D)

toc

D = DD;

tic

D = zKeepZeros(D,G);
nnz(D)

toc

% ------------------------------- Now test dropping zeros

D = DD;
G = zSparseRange(abs(E),1,2);           % cWW pairs

tic
D = D .* (G > 0);             
nnz(D)
toc

D = DD;
G = zSparseRange(abs(E),1,2,1);           % cWW pairs

tic
D = D .* G;             
nnz(D)
toc

D = DD;

tic
[i,j,d] = find(D);
[p,q] = find(G);
[c,k] = intersect([i j],[p q],'rows');
D = sparse(i(k),j(k),d(k));
nnz(D)
toc
