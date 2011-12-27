% pCompareSequencesFromFASTAFiles reads fasta files from loops extracted
% from 3D structures to identify groups that have the same sequences
% It ignores the flanking basepairs
% One should specify the loopType

disp('Make sure the Matlab current folder has a MotifLibrary in it');

SequenceSource = 1;               % 1 - Rfam, 2 - one Rfam, 3 - 3D structures
current = 1;

fs = 16;                          % font size for figures

% --------------------------------- Determine loop type

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  break
end

switch loopType,
case 'JL'
  R = 3;                      % for 3-way junctions, anyway
case 'IL'
  R = 2;                      % two rotations are computed
case 'HL'
  R = 1;
end

% --------------------------------- read sequence and model file names

switch SequenceSource
case 1,                          % many Rfam fasta files
  SeqPath = ['Sequences' filesep 'Rfam'];
  SeqFN = dir([SeqPath filesep '*.fasta']);
case 2,                          % one Rfam fasta file
  SeqFN(1).name = 'RF00431_IL_76_84_121_128.fasta';
  current = 1;
case 3,                          % fasta files from 3D structures
  clear SeqFN
  SeqPath = ['Sequences'];
  sequenceNameFile = [loopType '_Sequences.txt'];
  SeqNames = textread(['Sequences' filesep sequenceNameFile],'%s');
  for n = 1:length(SeqNames),
    SeqFN(n).name = SeqNames{n};
  end
end

NumGroups = length(SeqFN);

% diary pCompareSequencesFromStructures.txt
% diary pCompareSequencesFromRfam.txt

% ------------------------------- Loop through sequence files, get internal seq

clear Seq

for f = current:NumGroups,

  FASTA = zReadFASTA([SeqPath filesep SeqFN(f).name]);

  clear A

  for n = 1:length(FASTA),
    s = FASTA(n).Sequence;
    if R > 1,
      i = strfind(s,'*');
      A{1+(n-1)*R} = s([2:(i-2) i (i+2):(end-1)]);
      A{2+(n-1)*R} = s([(i+2):(end-1) i 2:(i-2)]);
    else
      A{n} = s;
    end
  end

  A = unique(A);

  for a = 1:length(A),
    Seq{f,a} = A{a};
  end

end

% -------------------------------- Loop through, compare sequences



[L,M] = size(Seq);

for n = 1:L,

if mod(n,50) == 0,
  fprintf('Working on fasta file %d\n', n);
end

 for m = (n+1):L,
   Diff{n,m} = [];
   for i = 1:M,
    if ~isempty(Seq{n,i}),
     for j = 1:M,
      if ~isempty(Seq{m,j}),
       if strcmp(Seq{n,i},Seq{m,j}),
        Diff{n,m} = [Diff{n,m} 0];
        a = strfind(Seq{n,i},'*');
        if (a > 2) && (a < length(Seq{n,i})-1),
         fprintf('Sequence %s is found in %s and %s\n', Seq{n,i}, SeqFN(n).name, SeqFN(m).name);
        end
       end
      end
     end
    end
   end
 end
end
