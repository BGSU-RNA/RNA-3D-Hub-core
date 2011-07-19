% zListNonredundantSet lists PDB ID's and chains that are meant to be non-redundant

Verbose = 1;

NRList = 'Nonredundant_4A_2010-05-19_list';
NRList = 'Nonredundant_3A_2010-05-19_list';
NRList = 'Nonredundant_1,5A_2010-05-19_list';

Filenames = zReadPDBList(NRList,1);

File = zAddNTData(NRList,0,[],1);

fid = fopen([pwd filesep 'Web' filesep 'AnalyzedStructures' filesep NRList '.txt'],'w');

Vers = num2str(File(1).ClassVersion);

fprintf(fid,'# PDB_ID_FR3D_Version_%s\n',Vers);

fprintf(fid,'PDB_ID\tChain(s)\n');

Chains = zBestChains(File,Verbose);

for f = 1:length(File),
  fprintf(fid,'%s\t',File(f).Filename);
  fprintf(fid,'%s',Chains{f});
  fprintf(fid,'\n');
end

fclose(fid);



















