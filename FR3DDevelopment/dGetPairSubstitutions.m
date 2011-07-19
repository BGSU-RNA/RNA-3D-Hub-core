
Verbose = 1;

f = 1;

%File = zAddNTData('Nonredundant_2009-05-14_list');
%File = zAddNTData({'2avy','2aw4'});

Filename = '2avy';
Filename = '2aw4';

[File,f] = zAddNTData(Filename,0,File,Verbose);   % load PDB data

[i,j] = find((File(f).Edge > 20) .* (File(f).Edge < 24));  % stacks
bp = 0;
goal = 'stacks';
[i,j] = find((File(f).Edge > 0) .* (File(f).Edge < 13));   % basepairs
bp = 1;
goal = 'basepairs';

if 0 > 1,
 clear SubsMatrices
 clear Counts
end

SubsMatrices{max(i),max(j)} = [];

length(i)

for k = 1:length(i),

  if isempty(SubsMatrices{i(k),j(k)}),

  fprintf('#%d Finding %s similar to %s%s-%s%s %s\n', k, goal, File(f).NT(i(k)).Base, File(f).NT(i(k)).Number, File(f).NT(j(k)).Base,  File(f).NT(j(k)).Number, zEdgeText(File(f).Edge(i(k),j(k))));

   clf
   VP.Sugar = 1;
   zDisplayNT(File(f),[i(k) j(k)],VP);
   drawnow


  [D,DGeom] = xPairSubstitutions(File,f,i(k),j(k));

  L = length(D(:,1));
  S = zeros(4,4);
  for c = 1:L,
    if D(c,3) < 1,
      S(D(c,1),D(c,2)) = S(D(c,1),D(c,2)) + 1 - D(c,3);
    end
  end
  S = S / sum(sum(S));                        % normalize

  fprintf('Found %d, substitution matrix is:\n', L);

  fprintf('Substitution matrix according to IsoDiscrepancy\n');
  100*S

  SubsMatrices{i(k),j(k)} = S;
  Counts(i(k),j(k)) = L;

  L = length(DGeom(:,1));
  S = zeros(4,4);
  for c = 1:L,
    if DGeom(c,3) < 0.5,
      S(DGeom(c,1),DGeom(c,2)) = S(DGeom(c,1),DGeom(c,2)) + 1 - 2*DGeom(c,3);
    end
  end
  S = S / sum(sum(S));                        % normalize

  fprintf('Substitution matrix according to Geometric Discrepancy\n');
  100*S

  if bp == 1,
    fprintf('pIsoScore\n');
    100*pIsoScore(File(f).Edge(i(k),j(k)),File(f).NT(i(k)).Code,File(f).NT(j(k)).Code,2)
  end

%  L = min(length(D(:,1)),500);
%  SubsData{i(k),j(k)} = D(1:L,:);
%  SubsData{j(k),i(k)} = D(:,[2 1 3]);

   if mod(k,100) == 0,
     save(['SubsMatrices' File(f).Filename goal],'SubsMatrices');
   end

if 0 > 1,
   clf
   subplot(2,1,1)
   hist(D(:,3))
   subplot(2,1,2)
   hist(DGeom(:,3))
   drawnow
   pause
end

  end
end

save(['SubsMatrices' File(f).Filename goal],'SubsMatrices');

break

c = [];
for f = 1:length(File),
  c = [c cat(2,File(f).NT.Code)];
end

for i = 1:4,
  count(i) = sum(c == i);
end

perc = count / sum(count);

perc = [0.2300    0.2647    0.3230    0.1823];

M = perc'*perc;

% divide SubsMatrices by M entry by entry to scale each entry according
% to how frequent the bases are in the 3D structures.  This will be some
% kind of "relative frequency".
