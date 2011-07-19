
load PDBInfo

for i = 1131:length(t(:,1)),
  File = zAddNTData(t{i,1},0,[],1);
%  File = zBackboneConformation(File,1);
%  zSaveNTData(File);
end
