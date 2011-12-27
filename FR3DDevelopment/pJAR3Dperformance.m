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
  R = 3;                      % for 3-way junctions, anyway
case 'IL'
  R = 2;                      % two rotations are computed
case 'HL'
  R = 1;
end

% --------------------------------- read sequence and model file names

sequenceNameFile = [loopType '_Sequences.txt'];
modelNameFile = [loopType '_Models.txt'];

SeqNames = textread(['Sequences' filesep sequenceNameFile],'%s');
ModNames = textread(['Models' filesep modelNameFile],'%s');

NumGroups = length(SeqNames);
NumModels = length(ModNames);

% ------------------------------- Count basepairs in each model
% ------------------------------- This SHOULD be done on the signature
% ------------------------------- saved in a data file - need to update

for g = 1:NumModels,
  m = ModNames{g};
  if isempty(strfind(m,'Helix')),
    m = strrep(m,'.txt','');
    i = strfind(m,'cWW');
    m(i(1):(i(1)+2)) = '---';
    m(i(end):(i(end)+2)) = '---';
    i = find( (m=='t') + (m=='c') );
    NumNoncWWBasepairs(g) = length(i);
  else
    NumNoncWWBasepairs(g) = 0;
  end
   fprintf('%3d %s\n', NumNoncWWBasepairs(g), ModNames{g});
end

% ------------------------------- Concatenate sequence files

% --------- Note: OwnGroup is defined below for the case in which
%           the groups of sequences correspond exactly to the group of
%           models.  When we run batches of Rfam sequences against models,
%           this will give the wrong idea; we'll need something else.

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

% ------------------------------- Rank of own model using MLPS and Percentiles

for n = 1:NumSequences,
  g = OwnGroup(n);
  mlps = MLPS(n,:,:);
  m = max(mlps,[],3);
  MLPSRank(n) = length(find(m > MLPS(n,g,1)));

  perc = Percentile(n,:,:);
  m = max(perc,[],3);
  PercRank(n) = length(find(m > Percentile(n,g,1)));
end

figure(1)
clf
hist(MLPSRank,30)
title('Rank of own group using max log prob score','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Rank_Own_Group_By_MLPS.png']);

figure(2)
clf
hist(PercRank,30)
title('Rank of own group using percentile score','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Rank_Own_Group_By_Percentile.png']);

figure(3)
clf
plot(MLPSRank+rand(1,NumSequences),PercRank+rand(1,NumSequences),'.');
xlabel('Rank of own model by max log prob score','fontsize',fs);
ylabel('Rank of own model by percentile','fontsize',fs);
m = max(max(PercRank),max(MLPSRank));
axis([0 m 0 m]);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Rank_By_Percentile_Against_Rank_By_MLPS.png']);

% ------------------------------------------- Focus on models with non-cWW


j = find(NumNoncWWBasepairs(OwnGroup) > 0); % sequences from models w/
                                            % more than flanking pairs

mlpsrank = MLPSRank(j);
percrank = PercRank(j);

figure(4)
clf
maxrank = 10;
n = hist(min(mlpsrank,maxrank),1:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(1:maxrank,n/sum(n));
axis([0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
xlabel('Rank of own model by max log prob score','fontsize',fs);
ylabel('Frequency','fontsize',fs);
title([num2str(length(mlpsrank)) ' individual sequences versus group models'],'fontsize',fs);
set(gca,'XTick',1:10)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','>=10'})
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Rank_By_MLPS_for_structured_sequences.png']);

figure(5)
clf
maxrank = 10;
n = hist(min(percrank,maxrank),1:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(1:maxrank,n/sum(n));
axis([0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
xlabel('Rank of own model by percentile score','fontsize',fs);
ylabel('Frequency','fontsize',fs);
title([num2str(length(percrank)) ' individual sequences versus group models'],'fontsize',fs);
set(gca,'XTick',1:10)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','>=10'})
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Rank_By_Percentile_for_structured_sequences.png']);



axis
plot(1:NumGroups, gmin, 'r.');
hold on
plot(1:NumGroups, gmed, 'b.');
plot(1:NumGroups, gmax, 'g.');
title('Min, median, max of percentiles of sequences against their own group','fontsize',fs);
xlabel('Group number','fontsize',fs)
ylabel('Min (red) median (blue) max (green)','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[loopType '_Own_Percentile_Min_Max_Med_By_Group.png']);

% -----------------------------------------------------------------------

