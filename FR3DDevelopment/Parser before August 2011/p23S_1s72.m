% p23S_1s72 loads the 23S structure and makes a model for the 23S

if exist('File'),
  [File,FIndex] = zAddNTData('1s72',2,File);
else
  [File,FIndex] = zAddNTData('1s72',2);
end

nMin = min(find(cat(1,File(FIndex).NT.Chain) == '0'));
nMax = max(find(cat(1,File(FIndex).NT.Chain) == '0'));

Mask1s7223S = zSimpleSecondaryStructure(File(FIndex)); % remove pseudoknots and most tertiary

F = File(FIndex);
F.Edge = File(FIndex).Edge .* (Mask1s7223S > 0);

zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax);
Node = pShiftNodeIndices(Node,nMin);

S = pTheoreticalAlignment(Node,1);

M = strrep(S{2},'-','');
M = strrep(M,'(','');
M = strrep(M,')','');
M = strrep(M,'<','');
M = strrep(M,'>','');
M = strrep(M,'{','');
M = strrep(M,'}','');

fprintf('Diagnostic:  Compare bases from structure (first line) with bases in experimental alignment:\n');
fprintf('%s seq: %s\n', File(FIndex).Filename, cat(2,File(FIndex).NT(nMin:nMax).Base));
fprintf('pTheoret: %s\n', M);
fprintf('\n');

fprintf('Diagnostic:  Compare experimental alignment with Java alignment.\n');
fprintf('Experimental alignment for %s:\n', File(FIndex).Filename);
fprintf('%s\n', S{1});
fprintf('%s\n', S{2});
fprintf('\n');

pWriteJavaNodeFile(File(FIndex),Node,4);

