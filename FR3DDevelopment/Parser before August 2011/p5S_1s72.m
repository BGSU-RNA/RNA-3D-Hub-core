% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if exist('File'),
  [File,FIndex] = zAddNTData('1s72',2,File);
else
  [File,FIndex] = zAddNTData('1s72',2);
end

nMin = zIndexLookup(File(FIndex),'1(9)');    % start of 5S chain
nMax = zIndexLookup(File(FIndex),'122(9)');  % end of 5S chain

i = nMin:nMax;

% focus on the 5S, rather than the whole 30S:

F      = File(FIndex);
F.NT   = F.NT(i);
F.Edge = F.Edge(i,i);
F.NumNT= length(i);
F.Distance = F.Distance(i,i);
F.BasePhosphate = F.BasePhosphate(i,i);

nMin = 1;
nMax = length(i);

Mask = zSimpleSecondaryStructure(F);         % remove pseudoknots and tert

F.Edge = F.Edge .* (Mask > 0);

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
fprintf('%s seq: %s\n', File(FIndex).Filename, cat(2,F.NT(nMin:nMax).Base));
fprintf('pTheoret: %s\n', M);
fprintf('\n');

fprintf('Diagnostic:  Compare experimental alignment with Java alignment.\n');
fprintf('Experimental alignment for %s:\n', File(FIndex).Filename);
fprintf('%s\n', S{1});
fprintf('%s\n', S{2});
fprintf('\n');

pWriteJavaNodeFile(File(FIndex),Node,4);
