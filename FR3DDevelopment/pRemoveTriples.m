% pRemoveTriples(File) simplifies base triples by keeping only the highest-priority non-crossing interaction; if a nucleotide makes all crossing interactions, the one with the lowest crossing number is kept

function [File] = pRemoveTriples(File, Verbose)

for f = 1:length(File),
  E = File(f).Edge;
  B = (abs(E) > 0) .* (abs(E) < 13);
  S = sum(B);                                 % number of basepairs
  t = find(S > 1);                            % NTs with more than 1
  for i = 1:length(t),
    j = find(B(t(i),:) > 0);                  % basepairs
    
    e = B(t(i),j);                            % interactions
    e = e + 100*(File.Crossing(t(i),j) > 0);  % penalize non-local ones
    [y,k] = min(e);                        % highest priority pair

    if Verbose > 0,
      fprintf('\nKeeping  %s%4s - %s%4s %4s Crossing %d\n', File(f).NT(t(i)).Base, File(f).NT(t(i)).Number, File(f).NT(j(k)).Base, File(f).NT(j(k)).Number, zEdgeText(File(f).Edge(t(i),j(k))),File(f).Crossing(t(i),j(k)));
    end

    for a = 1:length(j),
      if a ~= k,
        if Verbose > 0,
          fprintf('Removing %s%4s - %s%4s %4s Crossing %d\n', File(f).NT(t(i)).Base, File(f).NT(t(i)).Number, File(f).NT(j(a)).Base, File(f).NT(j(a)).Number, zEdgeText(File(f).Edge(t(i),j(a))),File(f).Crossing(t(i),j(a)));
        end

        File(f).Edge(t(i),j(a)) = 0;
        File(f).Edge(j(a),t(i)) = 0;
      end
    end
  end
end