MF = ['AltNodeTest2-NoE.txt'];
SeqFile = 'RF00001_reduced_labeled.txt';
NumSeq = 125;

Verbose = 0;
JAR3D_path                             % tell where JAR3D class files are 

PD = JAR3DMatlab.Align2(pwd,SeqFile,MF,NumSeq,0,15);
A = PD.sdata;
probsM = double(PD.probsM(1:size(PD.probsM),1:size(PD.probsM(1))));
Alig = Alignment.getAlignment(A,NumSeq);

pDisplayAlignment(Alig,probsM);