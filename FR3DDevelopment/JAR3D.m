% JAR3D.m runs the Java files and returns results to Matlab

% javaaddpath 'C:\Documents and Settings\zirbel\My Documents\JAR3D';
% javaclasspath              % show the Java search path within Matlab

% clear java                 % re-load Java classes

function [void] = JAR3D(FastaFile,ModelFile,NumSequences,DNA,Range)

if nargin < 3,
  NumSequences = 10;
end

if nargin < 4,
  DNA = 0;
end

if nargin < 5,
  Range = 20;
end

JAR3D_path                                % tell where JAR3D class files are 

C = ['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" ' FastaFile ' ' ModelFile ' ' num2str(NumSequences) ' ' num2str(DNA) ' ' num2str(Range)]

A = JAR3DMatlab.Align(pwd,'16S_sequence_from_2avy_1j5e.fasta','16S_from_2AVY.txt',2,0,15);


DNA = 0;
FASTA = [];

System.setProperty('user.dir','C:\Documents and Settings\zirbel\My Documents\JAR3D\bin');

JAR3D('16S_sequence_from_2avy_1j5e.fasta','16S_from_2AVY_with_motifs.txt',2,0,20);

sequenceData = Alignment.loadFastaColumnsDNA('C:\Documents and Settings\zirbel\My Documents\JAR3D\',0,0,DNA); 

sequenceData = Alignment.loadFastaColumnsDNA('C:\Documents and Settings\zirbel\My Documents\JAR3D\',0,0,DNA); 

			numSequences = (int)(Double.parseDouble(args[3]));
			range        = (int)(Double.parseDouble(args[5]));
			sequenceData = Alignment.doParse(sequenceData,numSequences,args[2],range);
			Alignment.displayAlignmentFASTA(sequenceData,numSequences);


system(C);







%['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" ' FastaFile ' ' ModelFile]);



% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequences_from_1j5e_2AVY.fasta 16S_from_2AVY.txt']);

% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequence_from_2avy.fasta 16S_from_2AVY.txt']);


