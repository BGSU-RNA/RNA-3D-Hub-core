% pHistogramMotifScores plots the logarithmic histogram of scores of many sequences against the model in Search, the scores of the actual sequences from the structures, and the scores of given sequences which may be a good match to the model

% Search is the .mat file for the motif
% aligcodes is the user-input set of sequences

% function [void] = pHistogramMotifScores(MotifName)

disp('Make sure the Matlab current folder is My Dropbox\BGSURNA\Motifs');

RNALett = 'ACGU';

Filenames = dir(['MotifLibrary' filesep 'Gro*']);

JAR3D_path;  loopType = 'IL'; 

for m = 2:length(Filenames),
  FN = ['MotifLibrary' filesep Filenames(m).name];
  load(FN,'Search','-mat')                             % Load search data

  fprintf('Working on %s\n', FN);

  Cand = Search.Candidates;
  [L,N] = size(Cand);                      % number of candidates
  N = N - 1;                               % number of nucleotides

  F.Edge = pConsensusInteractions(Search);     % find consensus interact list
  index = find(diag(fix(abs(F.Edge)),1)==1);
  Truncate = index+1;
  Search.Truncate = Truncate;
  T = Search.Truncate;

  MaxLoop = 8;                             % maximum number of internal IL NTs

  Skip = 0;

  MFN = Search.modelfilename;

  if ~isfield(Search,'logbins'),
    if N-4 > MaxLoop,                                % do random, not all
      FFN = ['Sequences' filesep 'Random_IL_' num2str(T-3) 'x' num2str(N-1-T) '.fasta'];
      if ~exist(FFN),
        fprintf('Writing Fasta file for random %d x %d internal loops\n', T-3, N-1-T);
        makeRandomSequencesFasta(FFN,T-3,N-1-T,4^MaxLoop);
      end

    else
      FFN=['Sequences' filesep 'All_IL_' num2str(T-3) 'x' num2str(N-1-T) '.fasta'];
      if ~exist(FFN),
        fprintf('Writing Fasta file for all %d x %d internal loops\n', T-3, N-1-T);
        makeAllSequencesFasta(FFN,T-3,N-1-T);
      end
    end

    FFN = strrep(FFN,'Sequences\','');

try
  tic
    scores = JAR3DMatlab.MotifParseSingle(pwd,FFN,MFN);
  toc
catch
  disp('Something did not parse right ******************');
  Skip = 1;
end

    if Skip == 0,
      [logcounts,logbins] = zLogCounts(exp(scores));
      Search.logcounts = logcounts;
      Search.logbins   = logbins;
      save(FN,'Search','-mat');
    end
  end

  if Skip == 0,
MFN
Search.ownsequencefasta

    selfscores = JAR3DMatlab.MotifParseSingle(pwd,Search.ownsequencefasta,MFN);
    selfscores = exp(selfscores);

selfscores

    Text = xFASTACandidates(Search.File,Search,0);
    selfseqs = Text(2:2:end);

selfseqs

    aligcounts = ones(1,length(selfscores));
    MinAligStructDist = zeros(length(selfscores),1);

    figure(1)
    clf
    zPlotMotifScores(selfseqs,aligcounts,selfscores,MinAligStructDist,Search.logcounts,Search.logbins)

    T  = strrep(Search.modelfilename,'.txt','');
    TT = strrep(T,'_','-');
    title(TT)
    orient portrait
    drawnow
    saveas(gcf,['Models' filesep T '.png']);

%    disp('Press a key to continue')
%    pause

  end

end

break


% --------------------------------------------- Determine sequences in Search

for c = 1:L,
  for n = 1:N,
    f = Cand(c,N+1);
    CandCodes(c,n) = File(f).NT(Cand(c,n)).Code;
  end
end

[codes,counts] = zUnique(CandCodes);

seqs = RNALett(codes);

if Verbose > 0,
  fprintf('Sequences from structures are:\n');
  seqs
end

% --------------------------------------- find consensus interaction list

[Edge,BPh,BR,Search] = pConsensusInteractions(Search);

if Verbose > 0,
  fprintf('Consensus interactions\n');
  full(Edge)
  full(BPh)
  full(BR)
end

% ------------------------------------------ consult sequence alignments

fprintf('Incorporating sequence data\n');
[File,aligseqs,aligcodes,aligcounts] = zGetVariantsFromAlignments(File,Search,Verbose);

clear File

% -------------------------------------- seqs in structures might not be in alig

StructAligDist = zDistanceHamming(codes, aligcodes); % distance to seqs in structures

MinStructAligDist = min(StructAligDist,[],2);  % distances
i = find(MinStructAligDist > 0);               % seq in struct not in alig
aligcodes = [aligcodes; codes(i,:)];           % append these codes
aligcounts = [aligcounts; zeros(length(i),1)]; % add in zero counts

% --------------------------------------- separate strands by a star

if exist(strandend),
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
  SS = Search;                                % store for later
  if Verbose > 0,
    fprintf('Loading stacking geometries\n');
  end
  load 2010-05-30_10_07_51-s35_NR_4A.mat
  Stacks{1} = Search;
  load 2010-05-29_10_24_40-s33_NR_4A_exclude_redundant.mat
  Stacks{2} = Search;
  load 2010-05-29_10_39_07-s55_NR_4A.mat
  Stacks{3} = Search;
  Search = SS;                                % restore Search
end

% ------------------------------------------ score all sequences

[allscores,ScoringMethodName,subscores,SSNames] = zScoreMotifSequences(allcodes,Search,Edge,BPh,Method,Stacks,BR);

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
for k = 3:5,                               % 

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
    T = sprintf('Sequence %4d %s Scores %16.16f %16.16f %16.16f', a, RNALett(allcodes(i(a),:)), allscores(i(a),1), allscores(i(a),2), allscores(i(a),3));

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


