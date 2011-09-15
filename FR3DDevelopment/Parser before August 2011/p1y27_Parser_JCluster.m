% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

clear Node
Filename = '1y27';

if exist('File'),
  [File,FIndex] = zAddNTData(Filename,0,File);
else
  [File,FIndex] = zAddNTData(Filename,0);
end

F = File(FIndex);

figure(1)
clf
zDisplayNT(F);

nMin = 1;
nMax = length(F.NT);

fprintf('Secondary structure with original interactions\n');
zSecondaryStructure(F,nMin,nMax);

% add a close cWW pair
F = pModifyEdge(F, '26(X)', '44(X)', 'cWW');

% remove interaction at the hairpin in loop 2:
F = pModifyEdge(F, '62(X)', '63(X)','');

% tertiary interactions between ends of the loops:
F = pModifyEdge(F, '33(X)', '66(X)','');
F = pModifyEdge(F, '34(X)', '65(X)','');
F = pModifyEdge(F, '35(X)', '64(X)','');
F = pModifyEdge(F, '37(X)', '61(X)','');
F = pModifyEdge(F, '38(X)', '60(X)','');

% remove interactions at the junction for now:
F = pModifyEdge(F, '46(X)', '53(X)','');
F = pModifyEdge(F, '49(X)', '76(X)','');
F = pModifyEdge(F, '20(X)', '76(X)','');
F = pModifyEdge(F, '21(X)', '75(X)','');
F = pModifyEdge(F, '22(X)', '52(X)','');
F = pModifyEdge(F, '23(X)', '46(X)','');

fprintf('Secondary structure with modified interactions\n');
zSecondaryStructure(F,nMin,nMax);

fprintf('Creating model nodes for this molecule\n');

Node = pMakeNodes(F,nMin,nMax);

Node(1).basehalfwidth = 20;

Node = pShiftNodeIndices(Node,nMin);

Node(17).subtype = 'XX';
Node(end).subtype = 'XXX';

Node(6).leftLengthDist = Node(5).leftLengthDist;
Node(6).rightLengthDist = Node(5).rightLengthDist;

Node(18).leftLengthDist = Node(8).leftLengthDist;

% Define junction cluster node ----------------------

n = 7;
   Node(n).type      = 'JunctionCluster';               % node type
   Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                           % state to state transitions
   Node(n).PIns	        = [0.05 0.95];     % when no previous state

   Node(n).Left(1,:)    = [1 2 3 4];       % nucleotides to use on left
   Node(n).LIP          = [1];             % probs for insertion possibs

   Node(n).Middle(1,:)  = [1 4 7 8];       % nucleotides to use in middle
   Node(n).MIP          = [1];             % probs for insertion possibs

   Node(n).Right(1,:)   = [1 2];           % nucleotides to use on right
   Node(n).RIP          = [1];             % probs for insertion possibs

   Node(n).Letters      ='CAAGCAGCGC';

   Node(n).LeftIndex    = 6;
   Node(n).MiddleIndex  = 32;
   Node(n).RightIndex   = 61;

   Node(n).IBases(1,:)  = [1 10];
   Node(n).Score(:,:,1) = pIsoScore(1,4,1);
   Node(n).IBases(2,:)  = [2 9];
   Node(n).Score(:,:,2) = pIsoScore(1,1,4);
   Node(n).IBases(3,:)  = [3 7];
   Node(n).Score(:,:,3) = pIsoScore(1,4,1);
   Node(n).IBases(4,:)  = [4 5];
   Node(n).Score(:,:,4) = pIsoScore(6,1,3);
   Node(n).IBases(5,:)  = [6 10]; 
   Node(n).Score(:,:,5) = pIsoScore(5,4,1);
   Node(n).IBases(6,:)  = [5 8];
   Node(n).Score(:,:,6) = pIsoScore(1,3,2);

   Node(n).Bl        = [10,15];
   Node(n).Bm        = [65,68];
   Node(n).Br        = [113];

% --------------------------------------------------------



Sequence = pReadFASTA('Purine_riboswitches.fasta',15,83);
Sequence = pSequenceConstraints(Node,Sequence);

Filename = 'Parser/Purine_riboswitches_junctioncluster';

disp('Parsing')
numtodo = 100;
k = 1;

for j=k:min(length(Sequence),(k+numtodo-1)),
  [Node,Sequence] = pParseSequences(Node,Sequence,j,j);
  k = j;
  Sequence = pDisplayMultipleAlignment(Node,Sequence);
  drawnow
  save(Filename,'Node','Sequence','k');
end

