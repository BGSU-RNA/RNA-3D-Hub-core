
Filenames = zReadPDBList('Nonredundant_4A_2010-05-19_list',1);
Filenames = zReadPDBList('Nonredundant_LT_4A_2011-01-07_list',1);

%mkdir('NonRedundantMatFiles');

for f = 1:length(Filenames),
  copyfile(['PrecomputedData\' Filenames{f} '.mat'],['NonRedundantMatFiles\' Filenames{f} '.mat'])
end
