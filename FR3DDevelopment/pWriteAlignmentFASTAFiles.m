% pWriteAlignmentFASTAFiles reads the alignments associated with each motif in the library and organizes them into a single fasta file

disp('Make sure the Matlab current folder is My Dropbox\BGSURNA\Motifs');

loopType = 'IL';                          % internal loops only

Param = [0 4];                            % Verbose, pair substitution method

Param = [1 4 0 1 10 1];                   % Verbose, etc.  6th is "use near"

% ----------------------------------------- Set variables

Verbose = Param(1);
Focus = loopType(1);
Types = {'HL','IL','JL'};                 % types of models we have

% ----------------------------------------- Read file names from the library

Filenames = dir(['MotifLibrary' filesep 'Gro*']);

keep = [];                               % of all models, which to keep

switch Focus,

case 'H',
  typ = 1;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:4),'LIB7') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'H'),
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

case 'I',
  typ = 2;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:3),'Gro') || (strcmp(Filenames(m).name(1:3),'LIB') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'I')), 
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

case 'J',
  typ = 3;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:3),'LIB') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'J'), 
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

end

Filenames = Filenames(find(keep));

% ----------------------------------------- Load sequence alignments,
%                                           write new fasta files

if ~exist('File'),
  File = [];
end

for m = 1:length(Filenames),
  MN = Filenames(m).name;
  FN = ['MotifLibrary' filesep MN '.mat'];
  load(FN,'Search','-mat')                             % Load search data

  fprintf('\nCurrent model is : %s %s\n', MN, Search.Signature);

  if Verbose > 0,
    fprintf('pMakeModelsFromLibrary:  Reading sequences from alignments for %s\n', MN);
  end

  MNSig = [MN '_' Search.Signature];
  FNSig = ['MotifLibrary' filesep MNSig '.mat'];

  Filenames(m).namesig = MNSig;

  % --------------------------------------- Get sequences from structures

  File = xDisplaySequenceForCandidates(File,Search);

  

  % --------------------------------------- Write sequences in FASTA format


end
