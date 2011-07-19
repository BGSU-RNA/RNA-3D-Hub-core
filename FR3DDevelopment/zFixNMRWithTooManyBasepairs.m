
load PDBInfo

i = find(n(:,3) == 0);               % structures without basepairs

for ii = 1:length(i),
  if n(i(ii),1) > 0,                   % x-ray structure
    F = zAddNTData(t{i(ii),1},4,[],0);   % try loading again
    E = abs(triu(F.Edge));
    j = sum(sum( (E > 0) .* (E < 13)));  % number of basepairs
    if j > 0,
      fprintf('File %s now has %4d basepairs ******************\n', t{i(ii),1}, j);
    end
  end
end
