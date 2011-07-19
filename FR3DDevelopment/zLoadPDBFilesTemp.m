
Names = zReadPDBList('Nonredundant_4A_2010-05-19_list.pdb');

for n = 2:length(Names),
  FN = ['G:\PDBFiles\' Names{n} '.pdb1'];
  if exist(FN),
    LFN = [pwd '\PDBFiles\' Names{n} '.pdb1'];
    fprintf('%d %s %s\n', n, FN, LFN);

    copyfile(FN,LFN);

    F = zAddNTData(Names{n},4,[],1);

    F
    F.NT(1)

    pause

    delete(LFN);
  end
end
    
