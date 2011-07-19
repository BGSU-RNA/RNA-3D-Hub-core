
% first, use FR3D to do symbolic searches for triples
% then load them here and make sense of them

load '2009-10-11_18_53_40-pair_pair_triple_search'

File = zAddNTData(Search.Filenames);

D = [];

for c = 1:length(Search.Candidates),
  f = Search.Candidates(c,4);               % file number
  i = Search.Candidates(c,1:3);             % indices of nucleotides
  D = [D; [File(f).Edge(i(1),i(2)) File(f).Edge(i(2),i(3)) File(f).Edge(i(1),i(3)) File(f).NT(i(1)).Code File(f).NT(i(2)).Code File(f).NT(i(3)).Code]];
end

D = sortrows(D,[1 2 4 5 6]);
