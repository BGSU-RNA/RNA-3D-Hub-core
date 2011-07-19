
load PDBInfo

i = find(n(:,3) == 0);               % structures without basepairs

names = [];

% avoid 2WYY.pdb1, it's huge and takes forever to read and adds nothing

for ii = length(i):-1:1,
  if n(i(ii),1) > 0 && ~strcmp(t{i(ii),1},'2WYY'),       % x-ray structure
    F = zAddNTData(t{i(ii),1},4,[],0);   % try loading again
    E = abs(triu(F.Edge));
    j = full(sum(sum( (E > 0) .* (E < 13))));  % number of basepairs
    if j > 0,
      fprintf('File %s now has %4d basepairs **********************\n', t{i(ii),1}, j);

      names = [names; t{i(ii),1}];
    end
  end
end
