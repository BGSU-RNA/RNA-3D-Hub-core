
Filenames = zReadPDBList('Nonredundant_2010-05-19_list',1);

mkdir('NRData');

for f = 1:length(Filenames),
  copyfile(['PrecomputedData\' Filenames{f} '.mat'],['NonRedundantMatFiles\' Filenames{f} '.mat'])
end
