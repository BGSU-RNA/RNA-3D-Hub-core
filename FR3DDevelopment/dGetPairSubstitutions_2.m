%the directory must be FR3D
f = 1;

% File = zAddNTData('Nonredundant_2009-05-14_list');
 File = zAddNTData({'2avy','2aw4'});

%[File,f] = zAddNTData('2avy',0,File,Verbose);   % load PDB data

[i,j] = find((File(f).Edge > 20) .* (File(f).Edge < 24));

sdi = zeros(length(i),1);

clear SubsData

for k = 1:length(i),

  fprintf('Finding stacks similar to %s%s-%s%s %s\n', File(f).NT(i(k)).Base, File(f).NT(i(k)).Number, File(f).NT(j(k)).Base, File(f).NT(j(k)).Number, zEdgeText(File(f).Edge(i(k),j(k))));

  D = xPairSubstitutions(File,f,i(k),j(k));
  %SubsData{i(k),j(k)} = D;
  sdi(k) = create_sdi(D);	%stacking_descrepancy_index
%  SubsData{j(k),i(k)} = D(:,[2 1 3]);

  close all

end

