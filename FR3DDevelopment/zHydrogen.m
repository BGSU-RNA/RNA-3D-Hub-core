[a,x] = xlsread('H_bonding_Atoms_from_Isostericity_Table.xls')

fid = fopen('zCheckHydrogen.m','w');       % open for writing

for r = 9:length(x(:,1)),
  