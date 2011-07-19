
FN = 'NonRedundant_2008_02_21_list';

% File = zAddNTData(FN,0,[],1);

fid = fopen([FN '_basepairs.txt'],'w');

DataHeader2 = sprintf('PDB_ID\tInteraction\tNucleotide_1_Base\tNucleotide_1_PDB_Number\tNucleotide_1_Chain\tNucleotide_1_Sequence_Position\tNucleotide_2_Base\tNucleotide_2_PDB_Number\tNucleotide_2_Chain\tNucleotide_2_Sequence_Position');

fprintf(fid,'%s\n',DataHeader2);

c = 1;

for f = 1:length(File),
  
  Chain = cat(2,File(f).NT.Chain);
  U     = unique(Chain);
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    chainoffset(u) = min(i);                    % first index in this chain
  end

  E   = File(f).Edge;

  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);
    e = E(i,:);
    j = find( (abs(e) > 0) .* (abs(e) < 13) );
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      fprintf(fid,'%s',strrep(DText{c},' ',''));

      c = c + 1;

    end
  end
end

fclose(fid);

