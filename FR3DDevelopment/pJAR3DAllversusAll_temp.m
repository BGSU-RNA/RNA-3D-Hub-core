% pJAR3DAllversusAll runs JAR3D on all motif sequences against all motif
% models.  One should specify the loopType and the SequenceSource

disp('Make sure the Matlab current folder has a MotifLibrary in it');

fs = 16;                          % font size for figures

% --------------------------------- Determine loop type

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  break
end

switch loopType,
case 'JL'
  R = 3;                      % three rotations, for 3-way junctions
case 'IL'
  R = 2;                      % two rotations are computed
case 'HL'
  R = 1;                      % only one "rotation"
end

% --------------------------------- Read sequence and model file names

sequenceNameFile = [loopType '_Sequences.txt'];
modelNameFile = [loopType '_Models.txt'];

SeqNames = textread(['Sequences' filesep sequenceNameFile],'%s');
ModNames = textread(['Models' filesep modelNameFile],'%s');

NumGroups = length(SeqNames);
NumModels = length(ModNames);

% ------------------------------- Concatenate sequence files

%           Note: OwnGroup is defined below for the case in which
%           the groups of sequences correspond exactly to the group of
%           models.  It won't work for Rfam sequences.

clear FASTA

for n = 1:NumGroups,
  newFASTA = zReadFASTA(['Sequences' filesep SeqNames{n}]);
  m = length(newFASTA);
  g = n*ones(m,1);                % group that these sequences belong to
  GroupSize(n) = m;
  if n == 1,
    FASTA = newFASTA;
    OwnGroup = g;
  else
    FASTA = [FASTA newFASTA];
    OwnGroup = [OwnGroup; g];
  end
end

NumSequences = length(FASTA);

% --------------------- Histogram number of instances in each group

for g = 1:NumGroups,
  GroupSize(g) = length(find(OwnGroup == g));
end

figure(1)
clf
hist(GroupSize,30);
title('Number of instances in each group','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Num_Instances.png']);

figure(1)
clf
semilogy(1:NumGroups,GroupSize,'.');
title('Number of instances versus group number','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Num_Instances_By_Group.png']);

% --------------------- Write sequences to one file for each rotation

r = 1;                          % first rotation
AllSequencesFile{1} = [loopType '_All_Sequences_' num2str(r) '.fasta'];
fid = fopen(['Sequences' filesep AllSequencesFile{1}],'w');
for n = 1:length(FASTA),
  fprintf(fid,'>%s\n',FASTA(n).Header);
  fprintf(fid,'%s\n',FASTA(n).Sequence);
end  
fclose(fid);

if R > 1,                       % more than one rotation, for IL, JL
  FASTA_R = FASTA;
end

for r = 2:R,
  for n = 1:length(FASTA_R),
    a = FASTA_R(n).Sequence;
    b = FASTA_R(n).Aligned;

    % fprintf('%s and %s become ', a, b);

    i = strfind(a,'*');
    a = [a((i(1)+1):end) '*' a(1:(i(1)-1))];
    i = strfind(b,'*');
    b = [b((i(1)+1):end) '*' b(1:(i(1)-1))];

    % fprintf('%s and %s.\n', a, b);

    FASTA_R(n).Sequence = a;
    FASTA_R(n).Aligned = b;

  end

  AllSequencesFile{r} = [loopType '_All_Sequences_' num2str(r) '.fasta'];
  fid = fopen(['Sequences' filesep AllSequencesFile{r}],'w');
  for n = 1:length(FASTA_R),
    fprintf(fid,'>%s\n',FASTA_R(n).Header);
    fprintf(fid,'%s\n',FASTA_R(n).Sequence);
  end  
  fclose(fid);
end

% ----------------------- Parse all sequences against all models, all rotations

JAR3D_path;  

fprintf('Parsing sequences against %4d models\n', NumModels);
for m = 1:length(ModNames),
  for r = 1:R,
    S = JAR3DMatlab.MotifParseSingle(pwd,AllSequencesFile{r},ModNames{m});
    Q = webJAR3D.getQuantilesB(S, SeqNames{m}(4:6), loopType);

    MLPS(:,m,r) = S;               % max log probability score for each seq
    Percentile(:,m,r) = Q;           % percentile of this score
  end

  if mod(m,50) == 0,
    fprintf('Parsed against %4d models so far\n', m);
  end
end

% MLPS(a,m,r) is the score of sequence a against model m using rotation r
% Percentile(a,m,r) is the percentile of sequence a against model m using r

% ------------------------------- Retrieve leave-one-out scores

load outfilename                % I think we need a new name for that file!

c = 0;
d = 1;
clear LeaveOneOut
for i = 1:length(FASTA),
  G = SeqNames{OwnGroup(i)};
  G = G(1:6);
  a = pCoreSequence(FASTA(i).Sequence);
  OwnScore(i) = MLPS(i,OwnGroup(i),1);

  for n = 1:length(group),
    b = pCoreSequence(sequence{n}{1});
    if strcmp(a, b) && strcmp(G, group{n}),
if scores(n) > OwnScore(i),
%      fprintf('%6s %20s %20s Core: %20s Scores: %10.6f  %10.6f\n', G, FASTA(i).Sequence, sequence{n}{1}, a, OwnScore(i), scores(n));
      c = c + 1;
end

if scores(n) < -900,
      fprintf('%6s %20s %20s Core: %20s Scores: %10.6f  %10.6f\n', G, FASTA(i).Sequence, sequence{n}{1}, a, OwnScore(i), scores(n));

end

      LeaveOneOut(i) = scores(n);
    end
  end

if length(LeaveOneOut) < i && OwnGroup(i) < 160,
%  fprintf('%6s %20s Core: %20s\n', G, FASTA(i).Sequence, a);

  OneVariant{d} = G(4:6);
  d = d + 1;
end

end

unique(OneVariant)


figure(3)
clf
plot(OwnScore(1:length(LeaveOneOut)),max(-30,LeaveOneOut),'.')
hold on
plot([-100 0],[-100 0],'r');
axis([-30 0 -30 0]);
xlabel('Score from own model');
ylabel('Score from leave-one-out model (min is -30)');

% ------ additional diagnostics concerning leave-one-out models

if 0 > 1,

for a = 1:length(FASTA),
  FASTASEQ{a} = FASTA(a).Sequence;
end
  
i = find(OwnGroup == 2);
unique(FASTASEQ(i))

c = 0;
clear JIMSEQ
for a = 1:length(group),
  if strcmp(group{a},'IL_002'),
    c = c + 1;   
    JIMSEQ{c} = sequence{a}{1};
  end
end
unique(JIMSEQ)

for a = 1:length(group),
  fprintf('%10s %20s Core: %20s\n', group{a}, sequence{a}{1}, pCoreSequence(sequence{a}{1}));
end


end

% ------------------------------- Histogram percentile against own model

for n = 1:NumSequences,
  g = OwnGroup(n);
  OwnPercentile(n) = Percentile(n,g,1);
end

figure(1)
clf
hist(OwnPercentile,30);
title('Percentile of each sequence against the model for its group','fontsize',fs);
ax = axis;
q = [0.999 0.99 0.98 0.97 0.96 0.95 0.9];
for m = 1:length(q),
  text(ax(1)+0.1*(ax(2)-ax(1)), (0.9-0.08*m)*ax(4), sprintf('%7.4f%% better than %7.3f', 100*sum(OwnPercentile > q(m))/NumSequences, q(m)));
end
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Own_Percentile_Histogram.png']);

figure(1)
clf
plot(OwnGroup,OwnPercentile,'.');
title('Own percentiles by group','fontsize',fs);
xlabel('Group number','fontsize',fs);
ylabel('Percentiles of sequences against their own group','fontsize',fs);

for g = 1:NumGroups,
  n = find(OwnGroup == g);        % sequences from group g
  gmed(g) = median(OwnPercentile(n));
  gmin(g) = min(OwnPercentile(n));
  gmax(g) = max(OwnPercentile(n));
end
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Own_Percentile_By_Group.png']);

figure(1)
clf
plot(1:NumGroups, gmin, 'r.');
hold on
plot(1:NumGroups, gmed, 'b.');
plot(1:NumGroups, gmax, 'g.');
title('Min, median, max of percentiles of sequences against their own group','fontsize',fs);
xlabel('Group number','fontsize',fs)
ylabel('Min (red) median (blue) max (green)','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Own_Percentile_Min_Max_Med_By_Group.png']);

% ------------------------------- Count models with better MLPS and Percentiles
% ------------------------------- Individual-group diagnostic

clear NumBetterScore
for n = 1:NumSequences,                     % loop through sequences
  g = OwnGroup(n);
  mlps = MLPS(n,:,:);                       % extract scores for this sequence
  m = max(mlps,[],3);

  NumBetterScore(n) = length(find(m > OwnScore(n)));
end

j = find(OwnGroup > 1);
% NumBetterScore = NumBetterScore(j);

figure(4)
maxrank = 10;
n = hist(min(NumBetterScore,maxrank),0:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(0:maxrank,n/sum(n));
axis([-0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
set(gca,'XTick',0:10)
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','>=10'})
title('Individual-group diagnostic');
xlabel(['Number of models that score higher than correct model']);
ylabel('Frequency');
























% ------------------------------- 

% -----------------------------------------------------------------------

break

% old code below

% ------------------------------- Run each file of sequences against all models

clear Seqs

for n = 1:length(SeqNames),
  Group{n}.Name = SeqNames{n};
  Group{n}.Number = n;

  clear NewSeq

SeqNames{n}

  FASTA = zReadFASTA(['Sequences' filesep SeqNames{n}]);
  
  NS = length(FASTA);


  if n == 1,
    Seqs = NewSeq;
  else
    Seqs = [Seqs NewSeq];
  end

end

% ------------------------------- Look up percentiles from mean scores in S

NumGroups = length(SeqNames);

if length(SeqNames) ~= length(S(:,1)),
  fprintf('Different number of sequence names and rows of parsing data!\n');
end

clear Q
for n = 1:NumGroups,
  s = [];
  for r = 1:R,
    s = [s; S(:,r+n*(r-1))];
  end
  FN = SeqNames{n}(4:6);
  sc = webJAR3D.getQuantilesB(s, FN, loopType);
  Q(:,n) = sc';
end

% ------------------------------- Store basic data about sequences

for n = 1:NumGroups,
  Sequence(n).Filename = SeqNames{n};
  Sequence(n).Group = n;
  Sequence(n).MeanMaxLogProbScore = reshape(S(n,:),R,NumGroups);
end

% ------------------------------- Look up percentiles for each MLP score

for n = 1:NumGroups,
  T = Sequence(n);
  for r = 1:R,
    for m = 1:NumGroups,
      Q = Sequence(m);
      sc = webJAR3D.getQuantilesB(T.MeanMaxLogProbScore(r,m), Q.Filename(4:6), loopType);
      Sequence(n).PercentileFromMean(r,m) = sc;
    end
  end

  Sequence(n).Percentile(1,n)

end

% ------------------------------- Histogram of MLPS for seqs against own model

ownMLPS = [];
ownPercentile = [];
for n = 1:NumGroups
  ownMLPS(n) = Sequence(n).MaxLogProbScore(1,n);
  ownPercentile(n) = Sequence(n).Percentile(1,n);
end
clf
hist(ownPercentile,30);

% ------------------------------- Save data

switch SequenceSource,
case 0,
  save(['pJAR3DAllversusAll_new_' loopType '.mat'],'Sequence');
case 1,
  save(['pJAR3DAllversusAll_new_' loopType '_Separated.mat'],'Sequence');
case 2,
  save(['pJAR3DAllversusAll_new_' loopType '_MSA.mat'],'Sequence');
end

