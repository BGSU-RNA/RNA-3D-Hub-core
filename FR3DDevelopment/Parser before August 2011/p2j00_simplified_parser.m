% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if exist('File'),
  [File,FIndex] = zAddNTData('2j00',2,File);
else
  [File,FIndex] = zAddNTData('2j00',2);
end

F = File(FIndex);

nMin = 1;
nMax = sum(cat(1,File(FIndex).NT.Chain)=='A');

Mask = zSimpleSecondaryStructure(File(FIndex)); % remove pseudoknots and most tertiary

F.Edge = File(FIndex).Edge .* (Mask > 0);

zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax);
Node = pShiftNodeIndices(Node,nMin);

pWriteJavaNodeFile(File,Node,4);

