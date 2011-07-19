% pMotifModel sets up nodes for a model of a motif

function [File,Node] = pMotifModel(File, Filename, StartIndex, EndIndex, Trunc)

if ~isempty('File'),
  [File,FIndex] = zAddNTData(Filename,2,File);
else
  [File,FIndex] = zAddNTData(Filename,2);
end

F = File(FIndex);

nMin     = zIndexLookup(F,StartIndex);     % 
nMax     = zIndexLookup(F,EndIndex);    %
Truncate = zIndexLookup(F,Trunc);     % where to cap with a **** hairpin

%fprintf('Original secondary structure\n');
%zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax,Truncate);

% get the indices down to reasonable numbers

a = Node(1).LeftIndex;
b = Node(1).RightIndex-length(Node);

for n = 1:length(Node),
  Node(n).LeftIndex  = Node(n).LeftIndex - a + 1;
  Node(n).RightIndex = Node(n).RightIndex - b + 1;
  Node(n).RightIndex = max(Node(n).RightIndex,max(Node(n).LeftIndex)+1);
  Node(n).MiddleIndex = Node(n).MiddleIndex - a + 1;
end

pWriteJavaNodeFile(File(FIndex),Node,5);

