

% 1s72, 2avy, 1j5e
% Produce a text output in N3 format of the FR3D-classified interactions in the given RNA molecule

% function [T] = zTextRNAOAssertions(File)

for f = 1:length(File),

clear T

r = 1;

RNAO_Base_Code{1} = 'RNAO_0000104';
RNAO_Base_Code{2} = 'RNAO_0000102';
RNAO_Base_Code{3} = 'RNAO_0000105';
RNAO_Base_Code{4} = 'RNAO_0000103';

% implement the following as a function, so it can accept negative arguments

T{r} = '@prefix rnao:     <http://purl.obofoundry.org/obo/rnao/>.';
r = r + 1;
T{r} = '@prefix rdfs:     <http://www.w3.org/2000/01/rdf-schema#>.';
r = r + 1;
T{r} = '@prefix ro:       <http://purl.obolibrary.org/obo/>.';
r = r + 1;
T{r} = '@prefix biotordf: <http://bio2rdf.org/>.';
r = r + 1;


T{r} = ['@prefix nt:       <http://bio2rdf.org/pdb:' File(f).Filename '_m1_c>.'];
r = r + 1;


T{r} = ['@prefix :         <http://rna.bgsu.edu/FR3D/rdf/' File(f).Filename '#>.'];
r = r + 1;
T{r} = '';
r = r + 1;

T{r} = [':m rdfs:label "This is a model of an RNA molecule based on experimental data from the pdb file' File(f).Filename ' ".'];
r = r + 1;
T{r} = [':m rdfs:label "Pdb descriptor: ' File(f).Info.Descriptor ' ".'];
r = r + 1;
T{r} = [':m :pdb_code "' File(f).Filename '".'];
r = r + 1;
T{r} = ':m rdfs:type rnao:RNAmolecule.';

% ---------------------------------------- introduce the nucleotides

for n = 1:length(File(f).NT),
  r = r + 1;
  T{r} = ['nt:' File(f).NT(n).Chain '_r' File(f).NT(n).Number ' rnao:sub_molecule_of :m.'];
end

% ---------------------------------------- label for each nucleotide

for n = 1:length(File(f).NT),
  r = r + 1;
%%%%%  T{r} = [':nt' num2str(n) ' rdfs:label "nucleotide ' num2str(n) ' in chain ' File(f).NT(n).Chain ' is ' File(f).NT(n).Base File(f).NT(n).Number '".'];
  T{r} = ['nt:' File(f).NT(n).Chain '_r' File(f).NT(n).Number ' rdfs:label "' File(f).NT(n).Base File(f).NT(n).Number '_r' File(f).NT(n).Chain '".'];
end

% ---------------------------------------- base identity of each nucleotide

for n = 1:length(File(f).NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rdf:type rnao:' RNAO_Base_Code{File(f).NT(n).Code} '.      # RNA base ' File(f).NT(n).Base];
end

% ---------------------------------------- nucleotide number from PDB file

for n = 1:length(File(f).NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rnao:nucleotide_number "' File(f).NT(n).Number '".'];
end

% ---------------------------------------- chain from PDB file

for n = 1:length(File(f).NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rnao:chain "' File(f).NT(n).Chain '".'];
end

% ---------------------------------------- covalent connectivity - pretend!

for n = 1:(length(File(f).NT)-1),
%  r = r + 1;
%  T{r} = [':nt ' num2str(n) 'rnao:three_prime_to_five_prime_to :nt' num2str(n+1) '.'];
end

% ---------------------------------------- pairwise interactions

E = abs(File(f).Edge);

[i,j] = find((E > 0) .* (E < 13));

[y,k] = sort(i);
i = i(k);
j = j(k);

for k = 1:length(i),
  r = r + 1;
  T{r} = ['nt:' File(f).NT(i(k)).Chain '_r' File(f).NT(i(k)).Number ' rnao:' zRNAOPairwiseInteractions(File(f).Edge(i(k),j(k))) ' nt:' File(f).NT(j(k)).Chain '_r' File(f).NT(j(k)).Number '.'];
end

% ---------------------------------------- stacking interactions

[i,j] = find((E > 20) .* (E < 24));

[y,k] = sort(i);
i = i(k);
j = j(k);

for k = 1:length(i),
  r = r + 1;
  T{r} = ['nt:' File(f).NT(i(k)).Chain '_r' File(f).NT(i(k)).Number ' rnao:' zRNAOPairwiseInteractions(File(f).Edge(i(k),j(k))) ' nt:' File(f).NT(j(k)).Chain '_r' File(f).NT(j(k)).Number '.'];
end

% ---------------------------------------- print the output

for r = 1:length(T),
  fprintf('%s\n', T{r});
end

% ---------------------------------------- write the output to a file

fid = fopen(['rdf' filesep File(f).Filename '_FR3D.n3'],'w');
for r = 1:length(T),
  fprintf(fid,'%s\n', T{r});
end
fclose(fid);

end

% zFTPToRNA([File(f).Filename '_FR3D.n3'],'/home/FR3D/AnalyzedStructures/1J1U');


