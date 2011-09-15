
File = zAddNTData('2avy');

clc
JAR3D_path                                % tell where JAR3D class files are 

% ----------------------------------------------- Fix up the file first

File(1) = pModifyEdge(File(1),'454','478',0);      % remove a cWW-cWW triple
File(1) = pModifyEdge(File(1),'113','353',-11);    % crystallographer fix
File(1) = pModifyEdge(File(1),'415','428',8);      % crystallographer fix

F = File(1);

% remove interactions which would necessitate a junction cluster

F = pModifyEdge(F,'959','984',0);
F = pModifyEdge(F,'959','1221',0);
F = pModifyEdge(F,'197','220',0);
F = pModifyEdge(F,'939','1375',0);

% remove an interaction on the left strand of a junction
% to restore this, allow for interactions within an initial node or
% allow a cluster to have zero length on one strand or another

F = pModifyEdge(F,'959','957',0);

if 10 > 1,
  F = File;
else
  F = xAnnotateWithKnownMotifs(File(1),1);
end

  F = xAnnotateWithKnownMotifs(F,1,0,{'2009-07-31_17_40_25-GU_packing_interaction_2avy.mat'});

Node = pMakeNodes(F,[1 4 2 1],'107','338');
Node = pShiftNodeIndices(Node,Node(1).LeftIndex-1);

p = (1./(1:120)).^2;
Node(1).leftLengthDist = p/sum(p);

Node(4).Delete = 0.104;
Node(5).Delete = 0.103;
Node(6).Delete = 0.102;
Node(7).Delete = 0.1;

ModelFile = ['16S_JAR3D_2avy_107-338.txt'];
pWriteJavaNodeFile(F,Node,4,ModelFile);
fprintf('Running JAR3D\n');


% JAR3D model needs to allow for more insertions; some sequences are much longer!

NumSequences = 100;

tt = cputime;

A = JAR3DMatlab.Align(pwd,'Knight-16S-soil-segments.fasta',ModelFile,NumSequences,2,15);

ttt = cputime;

fprintf('Alignment of %d sequences took %8.4f minutes, or %8.4f seconds per sequences\n', NumSequences, (ttt-tt)/60, (ttt-tt)/NumSequences);

Alig = Alignment.getAlignment(A,NumSequences);

fprintf('JAR3D alignment\n');
fprintf('%s\n', Alig.get(0));
for i = 1:NumSequences,
  fprintf('%s  %s\n', Alig.get(2*i), Alig.get(2*i-1));
end

fprintf('Alignment of %d sequences took %8.4f minutes, or %8.4f seconds per sequences\n', NumSequences, (ttt-tt)/60, (ttt-tt)/NumSequences);
