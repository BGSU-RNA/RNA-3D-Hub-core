% zDisplayManyAlignments makes one spreadsheet which displays Ryan's alignment, together with alignments from Jesse, JAR3D, ARTS, Needleman-Wunsch, SARSA, etc.

Molecule = '16S';

% read Jesse's alignment and put T.th first, E. coli second

[File2,i2,File1,i1] = zReadJesseAlignments('16S');

% I wasn't able to get this to work, besides, we need to use the new version:

rWriteAlignmentBPComparison(File1,1:length(File1.NT),File2,1:length(File2.NT),i1,i2);




% read in results of a JAR3D alignment of 1j5e and 2avy sequences, using 1j5e model

pGetAlignmentIndicesFromJAR3D

if 10 > 1,
  for m = 1:60,
    fprintf('JAR3D aligns %1s%4s with %1s%4s\n', File1.NT(i1(m)).Base, File1.NT(i1(m)).Number, File2.NT(i2(m)).Base, File2.NT(i2(m)).Number);
  end
end


S1 = cat(2,File1.NT.Base);
S2 = cat(2,File2.NT.Base);

[matches,i1,i2,s1,s2] = dNeedlemanWunsch(S1,S2,0.95,2);

if 10 > 1,
  for m = 1:60,
    fprintf('Needleman-Wunsch aligns %1s%4s with %1s%4s\n', File1.NT(i1(m)).Base, File1.NT(i1(m)).Number, File2.NT(i2(m)).Base, File2.NT(i2(m)).Number);
  end
end
