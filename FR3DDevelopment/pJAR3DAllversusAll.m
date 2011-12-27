% pJAR3DAllversusAll runs JAR3D on all motif sequences against all motif
% models.  One should specify the loopType and the SequenceSource

disp('Make sure the Matlab current folder has a MotifLibrary in it');

if ~exist('loopType'),
%  loopType = 'HL';
  loopType = 'IL';
end

if ~exist('SequenceSource'),
  SequenceSource = 1;                          % parse separate sequences?
end

% ----------------------------------------- Run JAR3D on these models

%clear java; JAR3D_path;  loopType = 'IL'; clc; S = JAR3DMatlab.MotifTest(pwd,loopType);
%clc; JAR3D_path;  
%JAR3DMatlab.MotifTest(pwd,loopType);

JAR3D_path;  

%S = JAR3DMatlab.MotifTest(pwd,loopType);

switch SequenceSource,
case 0,
  sequenceNameFile = [loopType '_Sequences.txt'];
case 1,
  sequenceNameFile = [loopType '_SeparatedSequences.txt'];
  disp('Make sure that you have run pSeparateSearches');
case 2,
  sequenceNameFile = [loopType '_SequencesMSA.txt'];
end

modelNameFile = [loopType '_Models.txt'];

S = JAR3DMatlab.MotifTestGeneral(pwd,loopType,sequenceNameFile,modelNameFile);

% entry S(i,2j-1) is the mean log likelihood score of sequence file i against model j
% entry S(i,2j)   is the mean log likelihood score of sequence file i against reversed model j

% [s,t] = size(S);
% switch loopType,
% case 'IL'
%   S = S(:,1:((2*t)/3));                 % remove extra columns
% case 'HL'
%   S = S(:,1:t);
% end

switch SequenceSource,
case 0,
  save(['pJAR3DAllversusAll_' loopType '.mat'],'S');
case 1,
  save(['pJAR3DAllversusAll_' loopType '_Separated.mat'],'S');
case 2,
  save(['pJAR3DAllversusAll_' loopType '_MSA.mat'],'S');
end

% ----------------------------------------- Display results graphically
break;
pDisplayModelScores
