
File = zAddNTData('1EC6');

nt = cat(1,File.NT.Center);
aa = cat(1,File.AA.Center);

D = zDistance(nt,aa);

[i,j] = find(D < 8);

File.NT(i(1)).Number
File.NT(i(1)).Base
File.NT(i(1)).Chain

% Now, do a FR3D search for this nucleotide and a neighboring one.
% Find the search save file in SearchSaveFiles
% in Matlab, load seachsavefilename

Search = rmfield(Search,'File');

xDisplayCandidates(File,Search)
