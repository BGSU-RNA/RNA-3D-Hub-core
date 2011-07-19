% shift indices from what you find in the .pdb file to correspond to
% letters found in the .fasta file

function [Node] = pShiftNodeIndices(Node,nMin)

for n=1:length(Node),
  Node(n).LeftIndex = Node(n).LeftIndex - nMin + 1;
  Node(n).RightIndex = Node(n).RightIndex - nMin + 1;
  if ~isempty(Node(n).MiddleIndex),
    Node(n).MiddleIndex = Node(n).MiddleIndex - nMin + 1;
  end
end
