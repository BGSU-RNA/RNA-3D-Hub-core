% needs to be changed to display nucleotide number instead of index, all the 
% way through


function [void] = pDisplayNodes(File,Node)

for n=1:length(Node),
  N = Node(n);

  fprintf('%3d %s',n, N.type);
  switch N.type,
  case 'Initial',
    fprintf(' Left param: %f Right param: %f\n', N.lpar, N.rpar);

  case 'Basepair',
    fprintf(' %c %4s %s %c %4s', N.LeftLetter, File.NT(N.LeftIndex(1)).Number, zEdgeText(N.Edge), N.RightLetter, File.NT(N.RightIndex(1)).Number);
    fprintf('  Left insertion %6.4f Right insertion %6.4f\n', mean(N.lpar(1:16)), mean(N.rpar(1:16)));

  case 'Hairpin'
    fprintf(' %s\n', N.subtype);

  case 'Motif'
    fprintf(' %d to %d on left, %d to %d on right\n',N.LeftIndex(1),N.LeftIndex(2),N.RightIndex(1),N.RightIndex(2));
    N.LeftLetter
    N.RightLetter
    Letters = [N.LeftLetter' N.RightLetter'];   % all letters interacting
    for i=1:length(N.IBases(:,1)),
      fprintf('     %c %d %s %c %d\n', Letters(N.IBases(i,1)), N.IBases(i,1), zEdgeText(N.Edge(i)), Letters(N.IBases(i,2)), N.IBases(i,2));
    end

  case 'Junction'
    fprintf('\n');

  case 'JunctionCluster'
    fprintf('\n    Left bases:   ');
    for i=1:length(N.LeftIndex(:)),
      fprintf(' %c %4s, ', File.NT(N.LeftIndex(i)).Base, File.NT(N.LeftIndex(i)).Number);
    end
    fprintf('\n');

    fprintf('    Middle bases: ');
    for i=1:length(N.MiddleIndex(:)),
      fprintf(' %c %4s, ', File.NT(N.MiddleIndex(i)).Base, File.NT(N.MiddleIndex(i)).Number);
    end
    fprintf('\n');

    fprintf('    Right bases:  ');
    for i=1:length(N.RightIndex(:)),
      fprintf(' %c %4s, ', File.NT(N.RightIndex(i)).Base, File.NT(N.RightIndex(i)).Number);
    end
    fprintf('\n');

  otherwise,
    fprintf('\n');
  end
end
