
File = zAddNTData('2aw4');
NTList1 = '1091  C 1100  G 1071';
NTList2 = 'A 2077  U 2243  U 2075';

Indices1 = zIndexLookup(File,NTList1);
Indices2 = zIndexLookup(File,NTList2);

[disc,SuperR,CC1,CC2] = xDiscrepancyForTriples(File,Indices1,File,Indices2);

figure(1)
clf

VP.Sugar = 1;
VP.SuperR = SuperR;
VP.CC1    = CC1;
VP.CC2    = CC2;
VP.Write  = 1;

zSuperimposeNucleotides(File,NTList1,File,NTList2,VP,3);

figure(2)
clf

[disc2,shift2,superR2] = zSuperimposeNucleotides(File,Indices1,File,Indices2);
