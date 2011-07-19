
Filenames = zReadPDBList('Nonredundant_2009-05-14_list');

mkdir('New_PDB_Folder');

for f = 1:length(Filenames),
  copyfile(['PDBFiles' filesep Filenames{f} '.pdb'],['New_PDB_Folder' filesep Filenames{f} '.pdb']);
end

for f = 1:length(Filenames),
  copyfile(['PDBFiles' filesep Filenames{f} '.pdb1'],['New_PDB_Folder' filesep Filenames{f} '.pdb1']);
end
