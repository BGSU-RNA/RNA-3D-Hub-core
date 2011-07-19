% 

if exist('File'),
  [File,FIndex] = zAddNTData('1s72',2,File);
else
  [File,FIndex] = zAddNTData('1s72',2);
end

F = File(FIndex);

nMin     = zIndexLookup(F,'1095');     % 
nMax     = zIndexLookup(F,'1261');     %
Truncate = zIndexLookup(F,'1100');     % where to cap with a **** hairpin

fprintf('Original secondary structure\n');
zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax,Truncate);

a = Node(1).LeftIndex;
b = Node(1).LeftIndex+10;              % 10 is the length of seqence

for n = 1:length(Node),
  Node(n).LeftIndex  = Node(n).LeftIndex - a + 1;
  Node(n).RightIndex = Node(n).RightIndex - b + 1;
  Node(n).MiddleIndex = Node(n).MiddleIndex - a + 1;
end

pWriteJavaNodeFile(File(FIndex),Node,5);

