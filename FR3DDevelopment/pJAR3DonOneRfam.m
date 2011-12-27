% pJAR3DonRfam runs JAR3D on loops extracted from Rfam.
% One should specify the loopType

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

SeqPath = ['Sequences' filesep 'Rfam'];
SeqFN = dir([SeqPath filesep '*.fasta']);

clear SeqFN
SeqFN(1).name = 'RF00431_IL_76_84_121_128.fasta';

modelNameFile = [loopType '_Models.txt'];
ModNames = textread(['Models' filesep modelNameFile],'%s');

NumGroups = length(SeqFN);
NumModels = length(ModNames);

% ------------------------------- Loop through sequence files, parse, tally

current = 1;

if current == 1,
  clear TopModel TopMLPS MedianPercentile
end

for f = current:NumGroups,

FASTA = zReadFASTA([SeqPath filesep SeqFN(f).name]);

% --------------------- Write sequences to one file for each rotation

r = 1;                          % first rotation
AllSequencesFile = [loopType '_All_Sequences_' num2str(r) '.fasta'];
fid = fopen(['Sequences' filesep AllSequencesFile],'w');
for n = 1:length(FASTA),
  fprintf(fid,'>%s\n',FASTA(n).Header);
  fprintf(fid,'%s\n',FASTA(n).Sequence);
end  

FirstSeq{f} = FASTA(1).Sequence;

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

  for n = 1:length(FASTA_R),
    fprintf(fid,'>%s\n',FASTA_R(n).Header);
    fprintf(fid,'%s\n',FASTA_R(n).Sequence);
  end  
end

fclose(fid);

% ----------------------- Parse all sequences against all models, all rotations

JAR3D_path;  

clear MLPS
clear Percentile

% fprintf('Parsing %d sequences against %4d models\n', R*length(FASTA), NumModels);
for m = 1:length(ModNames),
  try
    S = JAR3DMatlab.MotifParseSingle(pwd,AllSequencesFile,ModNames{m});
    Q = webJAR3D.getQuantilesB(S, ModNames{m}(4:6), loopType);
  catch
    S = -Inf*ones(R*NumModels,1);
    Q = zeros(R*NumModels,1);
  end

  for r = 1:R,
    a = 1+(r-1)*length(FASTA);
    b = length(FASTA)-1;
    MLPS(:,m,r) = S(a:(a+b));         % max log probability score for each seq
    Percentile(:,m,r) = Q(a:(a+b));   % percentile of this score
  end

  if mod(m,50) == 0,
%    fprintf('Parsed against %4d models so far\n', m);
  end
end

% -------------------------- Summarize diagnostics for parsing

% -------------------------- Record replications of each sequence

clear RepNum
for n = 1:length(FASTA),
  i = strfind(FASTA(n).Header,',');
  j = strfind(FASTA(n).Header,'instances');
  RepNum(1,n) = str2num(FASTA(n).Header((i+1):(j-1)));
end

% -------------------------- Maximum score over the rotations

MaxMLPS = MLPS(:,:,1);
MaxPercentile = Percentile(:,:,1);

for r = 2:R,
  MaxMLPS = max(MaxMLPS,MLPS(:,:,r));
  MaxPercentile = max(MaxPercentile,Percentile(:,:,r));
end

MaxMLPS = max(MaxMLPS,-1000);         % remove -Inf scores

WeightedAverageMLPS = (RepNum * MaxMLPS) / sum(RepNum);

[y,i] = sort(-WeightedAverageMLPS);

FirstPercentile{f} = MaxPercentile(1,i(1));     % 

for k = 1:20,                                   % store the top 20
  TopModel(f,k) = i(k);
  TopMLPS(f,k) = WeightedAverageMLPS(i(k));

  meds = [];
  for j = 1:length(RepNum),
    meds = [meds ones(1,RepNum(j))*MaxPercentile(j,i(k))];
  end
  MedianPercentile(f,k) = median(meds);
  NumRepetitions(f) = sum(RepNum);
end

fn = strrep(SeqFN(f).name,'.fasta','');
fn = fn(1:min(end,30));

mn = strrep(ModNames{i(1)}(1:min(end,40)),'.txt','');
mn = mn(1:min(end,30));

[mmp,w] = max(MedianPercentile(f,:));
mnt = strrep(ModNames{i(w)}(1:min(end,40)),'.txt','');
mnt = mnt(1:min(end,30));

fprintf('Group %4d of %4d, %-30s with %5d sequences, %20s, matches %-30s with median percentile %6.4f, first percentile %6.4f, but %-30s has median percentile %6.4f\n', f, NumGroups, fn, sum(RepNum), FirstSeq{f}(1:min(end,20)), mn, MedianPercentile(f,1), FirstPercentile{f}, mnt, mmp);

end

diary off

save('pJAR3DonRfam','TopModel','TopMLPS','NumRepetitions','FirstSeq','MedianPercentile','ModNames','SeqFN');

% ------------------------------------- Summarize the scores

figure(1)
clf
hist(MedianPercentile(:,1),30);
title('Median percentile of model with best MLPS','fontsize',fs);
set(gca,'fontsize',fs);
saveas(gcf,[loopType '_JAR3D_on_Rfam_best_MLPS.png']);

figure(2)
clf
[M,k] = max(MedianPercentile,[],2);
hist(M,30)
title('Highest median percentile','fontsize',fs);
set(gca,'fontsize',fs);
saveas(gcf,[loopType '_JAR3D_on_Rfam_best_percentile.png']);

p = 0.95;

i = find(MedianPercentile(:,1) > p);
[TopModels,NumTimes] = zUnique(TopModel(i,1));

for j = 1:length(TopModels),
  fprintf('Model %-50s had the best MLPS and scored above %7.2f%% %4d times\n', ModNames{TopModels(j)}, 100*p, NumTimes(j));
end

i = find(M > p);
clear TopMedianModel

for j = 1:NumGroups,
  [y,a] = max(MedianPercentile(j,:));
  TopMedianModel(j) = TopModel(j,a);
  TopMedian(j) = MedianPercentile(j,a);
end

p = 0.95;

k = find(TopMedian > p);

[TopModels,NumTimes] = zUnique(TopMedianModel(k)');

for j = 1:length(TopModels),
  fprintf('Model %-50s had the highest percentile and scored above %7.2f%% %4d times\n', ModNames{TopModels(j)}, 100*p, NumTimes(j));
end
  
