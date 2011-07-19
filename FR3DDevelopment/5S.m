
File = zAddNTData({'1s72','2aw4','2j01'},2);

%1s72 H. marismortui

nMin = zIndexLookup(File(1),'1(9)');    % start of 5S chain
nMax = zIndexLookup(File(1),'122(9)');  % end of 5S chain

[cat(2,File(1).NT(nMin:nMax).Base) ' >A H.m. from structure']

%2aw4 E. coli

nMin = zIndexLookup(File(2),'2(A)');    % start of 5S chain
nMax = zIndexLookup(File(2),'118(A)');  % end of 5S chain

[cat(2,File(2).NT(nMin:nMax).Base) ' >B E.coli from structure']

%2j01 T. thermophilus

nMin = zIndexLookup(File(3),'1(B)');    % start of 5S chain
nMax = zIndexLookup(File(3),'119(B)');  % end of 5S chain

[cat(2,File(3).NT(nMin:nMax).Base) ' >B T.th. from structure']

