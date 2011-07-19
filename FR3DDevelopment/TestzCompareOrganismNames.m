
load PDBInfo

k = find(n(:,4) > 1000);              % large structures only
t = t(k,:);

for kk = 1:length(k),
  fprintf('%s %s\n', t{kk,1}, t{kk,8});
end


L = length(t(:,1));

while 1 > 0,
  i = ceil(rand*L);
  j = ceil(rand*L);
  
  in = lower(t{i,8});
  jn = lower(t{j,8});

  CompareOrg = zCompareOrganismNames(in,jn);  % 1 if diff't, 0 if same

  fprintf('%s and %s are ', in, jn);

  if CompareOrg > 0,
    fprintf('different.\n');
  else
    fprintf('the same.\n');
  end
end
