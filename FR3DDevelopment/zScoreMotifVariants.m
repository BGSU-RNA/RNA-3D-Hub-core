% zScoreMotifVariants takes a set of candidates from a FR3D search, called the BuildSearch, generates all possible sequence variants, and scores them according to various pairwise scoring schemes, then displays the scores of all sequence variants from 3D structures and from corresponding sequence alignments

% Search can be a Matlab variable from xFR3DSearch or a filename
% It can also be an instance; first four letters are PDB ID, others are NTs
% Method is the method used in pIsoScore
% Display tells which type of display to use

% function [void] = zScoreMotifVariants(Search,Method,Display,Verbose)

if 0 > 1,                          % for when it is run as a function

if strcmp(class(Search),'char'),
  if ~isempty(strfind(lower(Search),'.mat')),
    load(Search)
  else
    [F,NTList,Indices] = zParseFilenameNucleotideString(Search);
    clear Search
    Search.Candidates = [Indices 1];
    Search.File       = F;
  end
end

if nargin < 2,
  Method = 2;
end

if nargin < 3,
  Display = 1;
end

if nargin < 4,
  Verbose = 1;
end

else

  Verbose = 1;
end

% ------------------------------------------- parameters for text labels

x = 10;
y = 10.^(-[10:0.5:24]);

for motif = 9:9,

if motif ~= 0,
  clear BuildSearch TestSearch UseAlignments
end

switch motif
case 0,                                      % don't load a new motif

case 1,
  load 2011-02-25_10_14_10-Sarcin_7_NR_NoRibo_LT_4A_2011-02-19.mat
  BuildSearch = Search;
  load 2011-02-25_18_02_03-Sarcin_7_NR_4A_2011-02-19.mat
  TestSearch = Search;
  MotifName = 'Sarcin 7-nucleotide';
  UseAlignments = {'2QBG','2AW7'};           % only attach to E. coli
  strandend = 4;                       % last nucleotide on first strand
case 2,
  load 2010-05-22_23_06_50-Sarcin_9_mixed_2aw4_2avy.mat
  MotifName = 'Sarcin 9-nucleotide';
case 3,
%  load 2010-05-22_23_38_17-C_loop_core_with_flanking_cWW_centroid_tuned_2aw4_2avy.mat
  load 2010-05-23_10_17_08-C_loop_core_with_flanking_cWW_2aw4_2avy_ACAUAU.mat
  load 2011-02-26_11_41_20-C_loop_core_with_flanking_cWW_NR_4A_ACAUAU.mat
  MotifName = 'C-loop';
case 4,
  load 2010-05-22_23_54_01-Kink_turn_highly_constrained.mat
  MotifName = 'Kink-turn';
case 5,
  load 2010-05-25_12_25_48-IL_tSH_tWH_tHS_2aw4_2avy
  load 2011-02-26_10_46_33-OCCBIO_motif_1S72_NR_4A.mat
  Cand = Search.Candidates;
  Cand = Cand(:,[2 3 4 7 8 9 11]);       % remove flanking cWW pairs, keep file
  Search.Candidates = Cand;
  MotifName = 'tSH-tWH-tHS IL';
  strandend = 3;                         % last nucleotide in first strand
case 6,
  load 2010-05-23_00_28_00-Helix_10_nucleotides_2avy_2aw4.mat
  MotifName = '8-nucleotide helix';
case 7,
  load 2010-08-13_19_33_15-Viroid_motif_focused.mat
  Cand = Search.Candidates;
  Cand = Cand(:,[2 3 4 7 8 9 11]);       % remove flanking cWW pairs, keep file
  Search.Candidates = Cand;
  MotifName = 'Viroid internal loop'
case 8,
  load 2010-07-27_16_33_56-Sarcin_7_mixed_NR_set_1S72.mat
  MotifName = 'Sarcin 7-nucleotide NR';
case 9,
  load 2009-08-03_16_42_13-Sarcin-Ricin_13_local_non-redundant.mat
  load 2011-02-26_12_38_20-Sarcin-Ricin_13_NR_4A.mat
  Cand = Search.Candidates;
  Cand = Cand(:,[2:6 9:12 end]);       % remove flanking cWW pairs, keep file
  Search.Candidates = Cand;
  MotifName = 'Sarcin 9-nucleotide';
  strandend = 5;                       % last nucleotide on first strand
  x = 5;
  y = 10.^(-[13:0.5:24]);
case 10,
  load 2010-11-07_11_06_27-cWW-tSH-tWH-tHS-cWW_NR4A_symbolic_tuned.mat
  strandend = 5;                       % last nucleotide on first strand
  MotifName = 'cWW-tSH-tWH-tHS-cWW';
  x = 10;
  y = 10.^(-[10:0.5:24]);
end

if ~exist('BuildSearch'),
  BuildSearch = Search;
  TestSearch  = Search;
end

clear Search

diary(['Motif scoring\' MotifName '.txt'])

fprintf('Working on %s\n', MotifName);

Method  = 2;                             % pIsoScore method to use
Verbose = 1;

RNALett = 'ACGU';

% --------------------------------------- find consensus interaction list

[Edge,BPh,BR,BuildSearch] = pConsensusInteractions(BuildSearch);

if Verbose > 0,
  fprintf('Consensus interactions\n');
  full(Edge)
  full(BPh)
  full(BR)
end

% ---------------------------------------- 

Cand = TestSearch.Candidates;

[L,N] = size(Cand);                      % number of candidates in TestSearch
N = N - 1;                               % number of nucleotides

ff = unique(Cand(:,end));                % file numbers needed for candidates

Filenames(ff) = TestSearch.Filenames(ff);    % file names needed

if ~exist('UseAlignments'),
  for f = 1:length(ff),
    UseAlignments{f} = Filenames{ff(f)};
  end
end

% --------------------------------------------- Load files and alignments

if ~exist('File'),                       % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);

% ---------------------------------------- Determine sequences in TestSearch

clear CandCodes

for c = 1:L,
  for n = 1:N,
    f = Cand(c,N+1);
    CandCodes(c,n) = File(f).NT(Cand(c,n)).Code;
  end
end

[codes,counts] = zUnique(CandCodes);

seqs = RNALett(codes);

if Verbose > 0,
  fprintf('Sequences from structures in TestSearch are:\n');
  seqs
end

% ------------------------------------------ consult sequence alignments

fprintf('Incorporating sequence data\n');
[File,aligseqs,aligcodes,aligcounts] = zGetVariantsFromAlignments(File,TestSearch,Verbose,UseAlignments);

clear File                                 % need to conserve memory

% -------------------------------------- seqs in structures might not be in alig

StructAligDist = zDistanceHamming(codes, aligcodes); % distance to seqs in structures

MinStructAligDist = min(StructAligDist,[],2);  % distances
i = find(MinStructAligDist > 0);               % seq in struct not in alig
aligseqs   = [aligseqs; RNALett(codes(i,:))];  % append these sequences
aligcodes  = [aligcodes; codes(i,:)];          % append these codes
aligcounts = [aligcounts zeros(1,length(i))];  % add in zero counts

% --------------------------------------- separate strands by a star

if exist('strandend'),
  clear star

  for i = 1:length(aligseqs(:,1)),
    star(i,1) = '*';
  end

  aligseqs = [aligseqs(:,1:strandend) star aligseqs(:,(strandend+1):end)];
end

% ------------------------------------------- generate all possible sequences

allcodes = zeros(4^N,N,'uint8');             % all sequences of length N
c = ones(1,N);                               % start with all 1's
allcodes(1,:) = c;                           % first row
r = 2;                                       % current row
a = N;                                       % start with last digit

while a > 0,
  if c(a) < 4,
    c(a) = c(a) + 1;                         % move to next sequence
    allcodes(r,:) = c;                       % store this one
    r = r + 1;
  else
    while a >= 1 && c(a) == 4,
      c(a) = 1;                              % roll this one over
      a = a - 1;
    end
    if a > 0,
      c(a) = c(a) + 1;  
      allcodes(r,:) = c;
      r = r + 1;
      a = N;
    end
  end
end

% ------------------------------------------ distances to aligned sequences

AligAllDist = zDistanceHamming(allcodes, aligcodes); % distance to seqs in structures

clear ii
clear jj

for k = 1:length(aligcodes(:,1)),
  [i,j] = find(AligAllDist(:,k) == 0);            % matched pairs

  ii(k) = i;                                      % which sequence in allcodes
                                                  % matches aligcode k
  jj(k) = k;
end

MinAligAllDist = min(AligAllDist,[],2);
clear AligAllDist

% --------------------------------------------- Load all stacks from NR dataset

if ~exist('Stacks'),
  if Verbose > 0,
    fprintf('Loading stacking geometries\n');
  end
  load 2010-05-30_10_07_51-s35_NR_4A.mat
  Stacks{1} = Search;
  load 2010-05-29_10_24_40-s33_NR_4A_exclude_redundant.mat
  Stacks{2} = Search;
  load 2010-05-29_10_39_07-s55_NR_4A.mat
  Stacks{3} = Search;
end

% ------------------------------------------ score all sequences

[allscores,ScoringMethodName,subscores,SSNames] = zScoreMotifSequences(allcodes,BuildSearch,Edge,BPh,Method,Stacks,BR);

clear Stacks

subscores = subscores(ii,:);               % we only need subscores for obs seq

save(['Motif_scoring_' MotifName], 'allscores', 'allcodes', 'ScoringMethodName', 'subscores');

% ------------------------------------------ find scores of observed seqs

fprintf('Finding distances between all codes and observed codes\n');

StructDist = zDistanceHamming(allcodes, codes); % distance to seqs in structures

MinAllStructDist = min(StructDist,[],2);     % min distance to *some* structure

% ----------------------------------------- Score by dist to structure too

allscores = [allscores(:,1:4) 10.^(-double(MinAllStructDist)-rand(size(MinAllStructDist)))];

[i,j] = find(StructDist == 0);            % matched pairs

for m = 1:length(i),
  score(j(m),:) = allscores(i(m),:);
end

% ------------------------------------------ find scores of alignment sequences

aligscore = zeros(length(aligcodes(:,1)),5);

for m = 1:length(ii),
  aligscore(m,:) = allscores(ii(m),:);
end

% ------------------------------------------ show subscores of those in seqs

[a,b] = sort(-aligcounts(jj));
ii = ii(b);
jj = jj(b);

figure(19)
clf

map = colormap;

nc = length(map(:,1));                        % number of colors

xx = 1:length(subscores(1,:));

for m = 1:length(ii),
  plot(xx,subscores(m,:),'linewidth',max(1,ceil(log(aligcounts(jj(m))))),'color',map(round(m*nc/length(ii)),:));
  hold on
end

title('Subscores of each sequence variant');

fprintf('Subscores for each observed sequence variant\n');

for a = 1:length(SSNames),
  fprintf('Subscore %3d is %s\n', a, SSNames{a});
end

fprintf('Subscores for sequences observed in alignments:\n');

for m = 1:length(ii),
  fprintf('%s ', aligseqs(jj(m),:));
  fprintf('%7.4f ', subscores(m,:));
  fprintf('\n');
end

clear subscores

% ------------------------------------------ distance from aligned to struct

AligDist = zDistanceHamming(aligcodes, codes); % distance to seqs in structures

MinAligStructDist = min(AligDist,[],2);     % min distance to *some* structure

% --------- Plot histograms of scores, dots for sequences from structures, aligs

%for k = 1:length(allscores(1,:)),          % loop through scoring methods
for k = 1:5,                               % 

  [logcounts,logbins] = zLogCounts(allscores(:,k),MinAllStructDist);

  figure(k)
  clf

  zPlotMotifScores(aligseqs,aligcounts,aligscore(:,k),MinAligStructDist,logcounts,logbins)


  % --------------------------------------- List scores to screen

  [yy,i] = sort(-allscores(:,k));          % sort by decreasing score

  M = 16384;                               % number of sequences to display
  M = length(allcodes(:,1));               % number of sequences to display
  M = 200;

  fprintf('Scoring method: %s\n', ScoringMethodName{k});

  g = 0;                                  % number of observed sequences found
  a = 1;                                  % current sequence
  while g < length(ii) && a < M,          % go until all observed are listed
    T = sprintf('Sequence %4d %s Scores %16.16f', a, RNALett(allcodes(i(a),:)), allscores(i(a),k));

    T = [T sprintf(' Distances %2d %2d', MinAllStructDist(i(a)), MinAligAllDist(i(a)))];

    if MinAligAllDist(i(a)) == 0,            % observed in alignment
      AADist = zDistanceHamming(allcodes(i(a),:), aligcodes); % distance to seqs in structures
      yy = find(AADist == 0);      
      g = g + 1;
      T = [T sprintf(' observed %5d times', aligcounts(yy))];
    end

    fprintf('%s\n', T);

    a = a + 1;
  end

  title(MotifName);
%  xlabel(['Motifs scored by ' ScoringMethodName{k}]);
  xlabel('Dots show number of times the variant was observed in an alignment');
  ylabel('Histogram shows scores of all sequences for this motif, colored by distance to a 3D sequence');

  saveas(gcf,['Motif scoring\' MotifName '_BPhGeom_Depth' num2str(k) '_BP' num2str(Method) '.pdf'],'pdf');
  saveas(gcf,['Motif scoring\' MotifName '_BPhGeom_Depth' num2str(k) '_BP' num2str(Method) '.png'],'png');

end

diary off


end
