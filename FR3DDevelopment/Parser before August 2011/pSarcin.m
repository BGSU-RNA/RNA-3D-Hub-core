% pSarcin sets up nodes for a model of the Sarcin internal loop

if exist('File'),
  [File,FIndex] = zAddNTData('1s72',2,File);
else
  [File,FIndex] = zAddNTData('1s72',2);
end

F = File(FIndex);

nMin     = zIndexLookup(F,'75(9)');     % 
nMax     = zIndexLookup(F,'106(9)');    %
Truncate = zIndexLookup(F,'82(9)');     % where to cap with a **** hairpin

fprintf('Original secondary structure\n');
zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax,Truncate);

a = Node(1).LeftIndex;
b = Node(1).LeftIndex+16;

for n = 1:length(Node),
  Node(n).LeftIndex  = Node(n).LeftIndex - a + 1;
  Node(n).RightIndex = Node(n).RightIndex - b + 1;
  Node(n).MiddleIndex = Node(n).MiddleIndex - a + 1;
end

pWriteJavaNodeFile(File(FIndex),Node,5);

