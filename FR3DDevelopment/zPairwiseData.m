
function [void] = zPairwiseData(File)

c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers

File.Distance = zMutualDistance(c,Inf); % compute distances < 16 Angstroms

[i,j] = find((File.Distance < 30) .* (File.Distance > 0));


for k = 1:length(i),
  fprintf('%s_%s|%s_%s|%0.4f|%d|%s\n', File.NT(i(k)).Number, File.NT(i(k)).Chain, File.NT(j(k)).Number, File.NT(j(k)).Chain, File.Distance(i(k),j(k)), abs(i(k)-j(k)), zEdgeText(File.Edge(i(k),j(k))));
end
