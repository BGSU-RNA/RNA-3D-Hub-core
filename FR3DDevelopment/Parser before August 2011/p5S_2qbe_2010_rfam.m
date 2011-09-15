% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

AlignmentList = 3;

MF = ['5S_2QBE_JAR3D'];
SeqFile = 'RF00001_seed.txt';
ModelStart = '2_A';                      % nucleotide 2, chain A
ModelEnd   = '118_A';
NumSeq = 712;
Verbose = 0;

if ~exist('File')
  File = zAddNTData('2qbe');
end

Verbose = 1;

JAR3D_path                             % tell where JAR3D class files are 

F = File(1);

Shift = zIndexLookup(File(1),ModelStart);


% ----------------------------------------------- JAR3D 1 alignment
%       NIH BISTI Method 1 scoring basepairs, no clusters, 
%       no extension of stems, no adjustment
%       of basepairing probabilities due to LR interactions

jj = 1;

if any(jj == AlignmentList),

  F = pModifyEdge(F,'30_9','52_9',0);       % remove nested triple
  F = pModifyEdge(F,'51_9','52_9',0);       % remove nested triple
  F = pModifyEdge(F,'36_9','47_9',0);
  F = pModifyEdge(F,'43_9','46_9',0);
  F = pModifyEdge(F,'78_9','79_9',0);


  JNode = pMakeNodes(F,[Verbose 4 0 0 0],ModelStart,ModelEnd);
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  Node1 = JNode;

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
  
  probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);

  F = File(1);

end

% ----------------------------------------------- JAR3D 2 alignment
%       NIH BISTI Method 1 scoring, clusters, no extension of stems, 
%       no adjustment
%       of basepairing probabilities due to LR interactions

jj = 2;

if any(jj == AlignmentList),
  JNode = pMakeNodes(F,[Verbose 4 0 0],ModelStart,ModelEnd); 
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  Node2 = JNode;

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
  
  probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);
end

% ----------------------------------------------- JAR3D 3 alignment
%        NIH BISTI Method 1 scoring with extension of stems, but no adjustment
%        of basepairing probabilities due to LR interactions

jj = 3;
%Insert dummy basepair to prevent extension of stem
F = pModifyEdge(F,'30_0','40_9',1);
i1 = zIndexLookup(F,'30_0');
i2 = zIndexLookup(F,'40_9');
F.Crossing(i1,i2)=10;

if any(jj == AlignmentList),
  JNode = pMakeNodes(F,[Verbose 4 1 0],ModelStart,ModelEnd);
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  Node3 = JNode;

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
  probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);

end

% ----------------------------------------------- JAR3D 4 alignment
%        NIH BISTI Method 1 scoring with extension of stems, and adjustment
%        of basepairing probabilities due to LR interactions and BPh

jj = 4;

if any(jj == AlignmentList),
  JNode = pMakeNodes(F,[Verbose 4 1 1],ModelStart,ModelEnd);
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);
end

% ----------------------------------------------- JAR3D 5 alignment
%      NIH BISTI Method 1 scoring with extension of stems, adjustment
%      of basepairing probabilities due to LR, BPh interactions, and GU packing

jj = 5;

if any(jj == AlignmentList),
  F = xAnnotateWithKnownMotifs(F,1,0,{'2009-07-31_17_40_25-GU_packing_interaction_2avy.mat'});

  JNode = pMakeNodes(F,[Verbose 4 1 1],ModelStart,ModelEnd);
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);
end

% ----------------------------------------------- JAR3D 6 alignment
%        NIH BISTI Method 1 scoring with extension of stems, adjustment
%        of basepairing probabilities due to LR, BPh interactions, GU packing,
%        and nothing but cWW and vague hairpins in extensible stems

jj = 6;

if any(jj == AlignmentList),
  F = xAnnotateWithKnownMotifs(F,1,0,{'2009-07-31_17_40_25-GU_packing_interaction_2avy.mat'});

  JNode = pMakeNodes(F,[Verbose 4 2 1],ModelStart,ModelEnd);
  T([1 jj+1],:) = pMakeNodesDiagnostics(F,JNode);
  JNode = pShiftNodeIndices(JNode,Shift);

  ModelFile = [MF num2str(jj) '.txt'];
  pWriteJavaNodeFile(File(1),JNode,4,ModelFile);
  PD = JAR3DMatlab.Align2(pwd,SeqFile,ModelFile,NumSeq,0,15);
  A = PD.sdata;
  probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
  Alig = Alignment.getAlignment(A,NumSeq);

  pDisplayAlignment(Alig,probsM,jj);
end




