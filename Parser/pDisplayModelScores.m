% pDisplayModelScores displays the results of pMakeModelsFromLibrary
% S(i,2k-1) tells the score that model k gave to sequences from model i
% S(i,2k)   tells the score that model k gave to reversed sequences from i
% We would hope that S(i,2i-1) would be the largest score in row i.

% simplealign.
% JAR3D_path
% SimpleAlign.editDist('ienarst','nensrt')
% SimpleAlign.editDist(ModelNames,ModelNames)

disp('Make sure the Matlab current folder has a MotifLibrary in it');

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  break
end

fprintf('Loop type is %s\n', loopType);

if ~exist('SequenceSource'),
  SequenceSource = 1;                          % parse separate sequences?
end

if ~exist('basepaironly'),
  basepaironly = 1;                        % only keep models with a basepair
end

% ---------------------------------------------- Load parsing data

if SequenceSource > 0,
  load(['pJAR3DAllversusAll_' loopType '_Separated'],'S');
else
  load(['pJAR3DAllversusAll_' loopType],'S');
end

% -------------------------------------------- Set up BaseURL for IL, HL

switch loopType,
case 'HL'
  BaseURL = 'http://rna.bgsu.edu/research/anton/share/Hairpins/';
case 'IL'
  if ~isempty(strfind(pwd,'ilmay')),
    BaseURL = 'http://rna.bgsu.edu/research/anton/share/ilmay6/';
  elseif ~isempty(strfind(pwd,'iljun')),
    BaseURL = 'http://rna.bgsu.edu/research/anton/share/iljun5/';
  else 
    BaseURL = 'http://rna.bgsu.edu/research/anton/share/iljun6/';
  end
otherwise
  BaseURL = '';
end

% ---------------------------------------- Read model names

n = 1;
fid = fopen(['Models' filesep loopType '_Models.txt'],'r');
if fid > 0,
  L = 1;
  while L > -1,
    L = fgetl(fid);                      % read a line
    if L > -1,
      Names{n} = L;
      switch loopType
      case 'IL'
        Names{n+1} = [L ' reversed'];    % sequences have been parsed twice
        n = n + 2;
      case 'HL'
        n = n + 1;
      end
    end
  end
  fclose(fid);
else
  fprintf('Could not open file of model names\n');
  break
end

switch loopType
case 'IL'
  t = length(Names);
  ModelNames = Names(1:2:(t-1));
case 'HL'
  ModelNames = Names;
end

% ---------------------------------------- Read sequence names

n = 1;

if SequenceSource > 0,
  fid = fopen(['Sequences' filesep loopType '_SeparatedSequences.txt'],'r');
else
  fid = fopen(['Sequences' filesep loopType '_Sequences.txt'],'r');
end

if SequenceSource < 3,
  clear SequenceNames

  if fid > 0,
    L = 1;
    while L > -1,
      L = fgetl(fid);                      % read a line
      if L > -1,
        SequenceNames{n} = L;
        n = n + 1;
      end
    end
    fclose(fid);
  else
    fprintf('Could not open file of sequence names, using model names instead\n');
    SequenceNames = ModelNames;              % if can't read sequence names file
  end
end

% ------------------------------------------ Read model signatures

try
  [modnum,NI,Sig,RSig] = textread(['Models' filesep loopType '_Signatures.txt'],'%s\t%s\t%s\t%s');
  for m = 1:length(NI),
    NumInst(m) = str2num(NI{m});
  end
catch
  fprintf('Signature file not found\n');
end

if length(Sig) ~= length(ModelNames),
  fprintf('pDisplayModelScores:  Different number of models and signatures\n');
end

% ---------------------------------- Remove helical models (needs improvement)

if SequenceSource == 1,
  ModelNames = ModelNames(1:(end-4));          % omit helices
end

% ---------------------------------- Remove models with no basepair beside
% ---------------------------------- flanking cWW

if basepaironly > 0,
  Keep = pMoreThanFlankingPairs(ModelNames,loopType);

  i = find(Keep);                  % indices of models to keep

  fprintf('Keeping %d models out of %d because they have more than just flanking cWW pairs\n', length(i), length(ModelNames));

  switch loopType
  case 'IL'
    j = sort([2*i-1 2*i]);           % indices of models and reversed models
    Names = Names(j);
    S = S(:,j);
  case 'HL'
    Names = Names(i);
    S = S(:,i);
  end

  ModelNames = ModelNames(i);
  Sig = Sig(i);
  RSig = RSig(i);
  NumInst = NumInst(i);

end

% --------------------------------------- Histogram number of instances

figure(66)
clf

maxinst = 40;
n = hist(min(NumInst(1:(end-4)),maxinst),1:maxinst);
bar(1:maxinst,n);
axis([0.5 maxinst+0.5 0 max(n)*1.1]);

xlabel('Number of instances');
ylabel('Frequency');
title('Number of instances in each group');
set(gca,'XTick',[1 2 3 4 5 6 7 8 10 15 20 25 30 35 40])
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','10','15','20','25','30','35','40+'})



% --------------------------------------- Find distances between groups

[p,GroupDist,GSize] = zOrderMotifs(ModelNames,Sig,RSig,1,0);

% --------------------------------------- List groups/models ordered by sim

if 10 > 1,
  for i = 1:length(ModelNames),
    j = p(i);

    fprintf('%4d\t%sGroup_%s.html\t%s\t%d\t\t%s\n', i, BaseURL, ModelNames{j}(4:6), ModelNames{j},NumInst(j),RSig{j});
                                        % paste this into google spreadsheet
  end
end

% --------------------------------------- List groups/models ordered by sim


if 10 > 1,
  fid = fopen('Motifs.html','w');

  fprintf(fid,'<html>\n');
  fprintf(fid,'<table border=1>\n');
  
  for i = 1:length(ModelNames),
    j = p(i);
%    j = per(i);                         % use order from pCompareTopQuantiles

    fprintf(fid,'<tr><td>%4d</td><td><a href="%sGroup_%s.html">Group %s</a></td><td><img src="Group_%s.png" width="150"></td><td>%s</td></tr>\n', i, BaseURL, ModelNames{j}(4:6), ModelNames{j}(4:6), ModelNames{j}(4:6), ModelNames{j});
                                        % paste this into google spreadsheet
    if 0>1 && i < length(ModelNames),
%    MinMax(per(i),per(i+1))
      if MinMax(per(i),per(i+1)) < 0.20,
        fprintf(fid,'<tr><td>Previous and next models have less than 0.20 overlap of top quantiles</td></tr>\n');
      end
    end

  end
  fprintf(fid,'</table>\n');
  fprintf(fid,'<html>\n');
  fclose(fid);
end

% ---------------------------------------- Extract model numbers

clear MN
clear GroupNum

for m = 1:length(ModelNames),
  tem = ModelNames{m};
  MN(m) = str2num(tem(4:6));              % model numbers as numbers
  GroupNum{m} = tem(4:6);                 % 3-digit group numbers as strings
  GNtoIndex(MN(m)) = m;
end

for j = 1:length(ModelNames),
  SLab{j} = ModelNames{j}(1:9);           % short labels
  Lab{j} = ModelNames{j};                 % full labels - just copies ModelNames
end

% ---------------------------------------- Remove sequences corresponding
% ---------------------------------------- to models with no extra basepairs

if basepaironly > 0,
  Keep = zeros(1,length(SequenceNames)); % initially, don't keep any
  for i = 1:length(SequenceNames),
    group = SequenceNames{i}(4:6);       % extract group number
    if ismember(group,GroupNum),
      Keep(i) = 1;
    end
  end

  i = find(Keep);

  fprintf('Keeping %d sequences out of %d because they have more than just flanking cWW pairs\n', length(i), length(SequenceNames));

  SequenceNames = SequenceNames(i);
  S = S(i,:);
end

% ---------------------------------------- Determine size of score matrix

[s,t] = size(S);                         % matrix of scores

% for IL, entry S(i,2j-1) is the mean log likelihood score of sequence file i against model j
% for IL, entry S(i,2j)   is the mean log likelihood score of sequence file i against reversed model j

% ------------------------- Recast scores as relative to max across all models
% ------------------------- Provide diagnostics on scores

RelSeq = [];
Confusion = zeros(t,t);                  % Confusion(i,j) 
MatchRank = zeros(1,s);                  % rank of own model, from 1 on down
GeoDisc   = [];                          % Geometric Discrep between model
                                         % and higher ranked models

counter  = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;

JAR3D_path;  

% for internal loops, entry RelSeq(i,2j-1) is how much below the the mean log likelihood score of sequence file i against model j is compared to the max in this row
% for internal loops, entry RelSeq(i,2j)   is how much below the the mean log likelihood score of sequence file i against reversed model j is compared to the max in this row

for i = 1:s,                             % loop through all sequence files
  [m,k] = max(S(i,:));                   % max score over all models
  RelSeq(i,:) = S(i,:) - m;              % subtract max score from this row
  a = find(RelSeq(i,:) == 0);            % maximizers of the score, may be >1
  [b,k] = sort(-RelSeq(i,:));            % sort models from best to worst

  if strcmp(loopType,'IL') && 2*s == t,  % IL have 2 models for each sequence
    if RelSeq(i,2*i-1) == 0,               % own model gives highest score
      counter = counter + 1;               
    end
    if RelSeq(i,2*i-1) >= -2,              % own model gives good score
      counter2 = counter2 + 1;
    end
    if RelSeq(i,2*i-1) < -5,               % own model gives bad score
      counter3 = counter3 + 1;
    end

    if 0 > 1 && RelSeq(i,2*i-1) < -0.0000001, % flag models with worse scores
      fprintf('Sequences from %s get scored better by another model\n',ModelNames{i});
      fprintf('%-50s self  score: %9.4f ',ModelNames{i},S(i,2*i-1));
      fprintf('%sGroup_%s.html\n', BaseURL, ModelNames{i}(4:6));
      for aa = 1:length(a),
        fprintf('%-50s gives score: %9.4f ',Names{a(aa)},S(i,a(aa)));
        fprintf('%sGroup_%s.html\n', BaseURL, Names{a(aa)}(4:6));
      end
      fprintf('\n');
    end
  end

  ownmodelbest = 0;                         % indicator whether
  for aa = 1:length(a),                     % loop through top models
    if strcmp(Names{a(aa)}(4:6),SequenceNames{i}(4:6)), % model in same group does best
      ownmodelbest = 1;
    end
  end

  if ownmodelbest > 0,
    counter4 = counter4 + 1;
    MatchRank(i) = 1;                       % own model gets rank 1
  end

  if ownmodelbest == 0,         % display when different
    
    seq = '';
    fid = fopen(['Sequences' filesep SequenceNames{i}],'r');
    if fid > 0,
      L = fgetl(fid);                      % read a line
      seq = fgetl(fid);                    % use the second line
      fclose(fid);
    end

    fprintf('%s sequence(s) %s get scored better by a different model\n',SequenceNames{i},seq);

    j = 1;
    confuse = zeros(1,t);
    Models = [];

    while ~strcmp(Names{k(j)}(4:6),SequenceNames{i}(4:6)) && j < length(k),
      if mod(k(j),2) == 0 && strcmp(loopType(1:2),'IL'),
        tem = [Names{k(j)}(1:7) RSig{k(j)/2} '(R)'];
        tem(10) = ' ';
      else
        tem = Names{k(j)};
        tem = tem(1:(end-4));
      end
      fprintf('%-50s gives score: %9.4f ',tem,S(i,k(j)));
      if exist(['Models' filesep 'Emp. Distributions' filesep Names{k(j)}(1:6) '.txt'])==2,
        sc = webJAR3D.getQuantilesB(S(i,k(j)), Names{k(j)}(4:6), loopType);
      else
        sc = NaN;
      end
      fprintf('and quantile %9.4f ', sc);

      fprintf('%sGroup_%s.html\n', BaseURL, Names{k(j)}(4:6));

      confuse(1,k(j)) = confuse(1,k(j)) + 1;

      Models = [Models GNtoIndex(str2num(Names{k(j)}(4:6)))];

      j = j + 1;
    end

    MatchRank(i) = j;                       % own model gets rank j

    w = GNtoIndex(str2num(Names{k(j)}(4:6)));
    
    for m = 1:length(Models),
      GeoDisc = [GeoDisc GroupDist(w,Models(m))];   % append discrepancies with
    end                                             % higher ranked models
    GeoDisc = [GeoDisc 0];                          % current model

    fprintf('%-50s gives score: %9.4f ',Names{k(j)}(1:(end-4)),S(i,k(j)));
    if exist(['Models' filesep 'Emp. Distributions' filesep Names{k(j)}(1:6) '.txt'])==2,
      sc = webJAR3D.getQuantilesB(S(i,k(j)), Names{k(j)}(4:6), loopType);
    else
      sc = NaN;
    end

    fprintf('and quantile %9.4f ', sc);
    fprintf('%sGroup_%s.html\n', BaseURL, SequenceNames{i}(4:6));
    fprintf('\n');

    Confusion(k(j),:) = Confusion(k(j),:) + confuse;  % record under this group

  end

end

if 2*s == t,
  fprintf('%3d out of %d sequences (%7.4f%%) are given the highest score by their own forward model\n', counter, s, 100*counter/s);
  fprintf('%3d out of %d sequences (%7.4f%%) are given the highest score by a model from their group\n', counter4, s, 100*counter4/s);
  fprintf('%3d out of %d sequences (%7.4f%%) are given a good score by their own forward model\n', counter2, s, 100*counter2/s);
  fprintf('%3d out of %d sequences (%7.4f%%) are given a very poor score by their own forward model\n', counter3, s, 100*counter3/s);
else
  fprintf('%3d out of %d sequences (%7.4f%%) are given the highest score by a model from their group\n', counter4, s, 100*counter4/s);
end

% ------------------------------------- Histogram ranks of matches

figure(88)
clf
maxrank = 10;
n = hist(min(MatchRank,maxrank),1:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(1:maxrank,n/sum(n));
axis([0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);

xlabel('Rank of own model');
ylabel('Frequency');
switch SequenceSource
case 0,
  title([num2str(s) ' group sequences versus group models']);
case 1,
  title([num2str(s) ' individual sequences versus group models']);
case 2,
  
end

set(gca,'XTick',1:10)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','>=10'})

if SequenceSource > 0,
  title(['Ranks of ' num2str(s) ' individual sequences across ' num2str(length(ModelNames)) ' ' loopType ' group models']);
else
  title(['Ranks of ' num2str(s) ' group sequence sets across ' num2str(length(ModelNames)) ' ' loopType ' group models']);
end

fprintf('%d out of %d (%7.4f%%) sequences have the correct model in the top 5\n', sum(n(1:5)), sum(n), 100*sum(n(1:5))/sum(n));

% ------------------------------------- Histogram discrepancy with higher
% ------------------------------------- ranked models

figure(99)
clf
maxrank = 10;
hist(min(GeoDisc,10),20);
%axis([0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
xlabel('Geometric discrepancy between model and equal or higher ranked models');
ylabel('Frequency');

if SequenceSource > 0,
  title(['Ranks of ' num2str(s) ' individual sequences across ' num2str(length(ModelNames)) ' ' loopType ' group models']);
else
  title(['Ranks of ' num2str(s) ' group sequence sets across ' num2str(length(ModelNames)) ' ' loopType ' group models']);
end

% ------------------------------------- List quantiles of sequences
% ------------------------------------- against their own models

clear OwnGroupScore

for i = 1:length(SequenceNames),
  group = SequenceNames{i}(4:6);
  k = find(ismember(GroupNum,group));

  seq = '';
  fid = fopen(['Sequences' filesep SequenceNames{i}],'r');
  if fid > 0,
    L = fgetl(fid);                      % read a line
    seq = fgetl(fid);                    % use the second line
    fclose(fid);
  end

  if isempty(k) || length(k) > 1,
    fprintf('%-50s has 0 or more than one model for its group\n');
  else

    if exist(['Models' filesep 'Emp. Distributions' filesep SequenceNames{i}(1:6) '.txt'])==2,
      OwnGroupScore(i) = webJAR3D.getQuantilesB(S(i,2*k-1), group, loopType);
    else
      OwnGroupScore(i) = NaN;
    end
  
%    fprintf('%-50s %-12s has quantile score %9.4f from the model for its group\n', SequenceNames{i}, seq, OwnGroupScore(i));
  end
end

fprintf('\n\n');

[y,k] = sort(OwnGroupScore);

for i = 1:length(SequenceNames),
  if OwnGroupScore(k(i)) < 0.8,
    fprintf('%-50s has quantile score %9.4f from the model for its group\n', SequenceNames{k(i)}, OwnGroupScore(k(i)));
  end
end

% ------------------------------------- Plot score confusion matrix
% ------------------------------------- Doesn't seem all that useful
% ------------------------------------- Only a few sets of sequences get
% ------------------------------------- scored better, and then by many models

if 10 > 1,

figure(77)
clf

D = 1./(1+Confusion+Confusion');
p = zOrderbySimilarity(D);

pcolor(Confusion(p,p));
hold on
axis ij
shading flat
caxis([0 10]);
colormap('default')
map = colormap;
map = map(8:56,:);
colorbar('eastoutside')

break
end


% ------------------------------------- Label lengths

maxlen = 0;
for i = 1:s,
  maxlen = max(maxlen,length(Lab{i}));
end
for i = 1:maxlen,
  blank(i) = ' ';
end
for i = 1:s,
  FL = blank;
  FL(1:length(Lab{i})) = Lab{i};
  FullLab{i} = FL(7:end);
end

% ------------------------------------ scores relative to max in column

RelModel = [];

for j = 1:t,
  RelModel(:,j) = S(:,j) - max(S(:,j));  % subtract the maximum; make it zero
end

% ------------------------------------ 

figure(10)
clf

p = 1:s;
q = 1:t;

T = RelSeq(p,q);
T = RelSeq;
T(s+1,t+1) = 0;
pcolor(T);
hold on
axis ij
shading flat
axis([1 2*s+1 1 s+1]);
caxis([-8 0]);
colormap('default')
map = colormap;
map = map(8:56,:);
%map = [1 1 1; map];
map = [map; 1 1 1];

q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];          % how to permute models, keep together
  plot([2*i 2*i+2 2*i+2 2*i 2*i]-1,[i i i+1 i+1 i],'k');
  hold on
end

colormap(map);
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',FullLab(p),'FontSize',5,'FontName','FixedWidth')
title('Scores of sequences against models and reversed models, relative to max score for those sequences');
xlabel('Models and reversed models, paired, in the order on the y axis');

% ------------------------------------ Distances between sequence files

D = [];
D(s,s) = 0;

for i = 1:s,
  for j = 1:s,
    ii = [2*i-1 2*i];                % indices of forward and reversed
    jj = [2*j-1 2*j];
    D(i,j) = min(min([max([0 0],(S(i,ii)-S(i,jj))) max([0 0],(S(j,jj)-S(j,ii)))]));
  end
end

% --------------------------------- Graph of distances between sequence files

figure(11)
clf

p = zOrderbySimilarity(D);

T = D(p,p);
T(s+1,s+1) = 0;
pcolor(T);
hold on
axis ij
shading flat
caxis([0 8]);
colormap('default')
colorbar('eastoutside');
title('Similarity of models according to how they score each others sequences');

% ========================================== organize instances by model
% ========================================== should work with separated
% ========================================== sequences and instances

for m = 1:length(ModelNames),
  tem = ModelNames{m};
  fprintf('%4d %s\n', m, tem);
end

NumMod = length(ModelNames);                % number of models

MD = [];                                    % distance between models

DD = abs(zDistance(max(-4,RelSeq)));     % how similarly sequences are scored

MNP = MN;                                % nominal order

for i = 1:NumMod,                        % loop through models
  k = find(MN == MN(i));                 % sequences corresponding to model i

  ip = zOrderbySimilarity(DD(k,k));      % instance perm w/in model
                                         % order instances in model by simil.

  DDD = abs(zDistance(S(k,2*k-1)));      % distance acc. to scores w/in model

  ip = zOrderbySimilarity(DDD);          % instance perm w/in model
                                         % order instances in model by simil.
if 0 > 1
  S(k,2*k-1)
  DDD

  figure(41)
  clf
  pcolor(DD(k(ip),k(ip)));
  caxis([0 20]);
  axis ij
  colorbar('eastoutside');
  fprintf('Similarity of scores within %4d instances in model %s\n', length(k), ModelNames{k(1)});
  drawnow
  pause
end

  MNP(1,k(ip)) = 1:length(k);                % reorder instances within model

if 0 > 1,
  figure(41)
  clf
  [y,ip] = sort(MNP(1,k));
  pcolor(DD(k(ip),k(ip)));
  pause
end

end

% -------------------------------------------- distance between models

for i = 1:NumMod,
  k = find(MN == MN(i));                 % sequences corresponding to model i
  for j = (i+1):NumMod,
    m = find(MN == MN(j));
    MD(i,j) = mean(mean(DD(k,m)));           % how similarly all seqs score
    MD(j,i) = MD(i,j);
  end
end

mp = zOrderbySimilarity(MD);                 % order in which to put models

fprintf('Models ordered by how they score sequences\n');
for m = 1:length(ModelNames),
  tem = ModelNames{mp(m)};
  fprintf('%4d %s\n', m, tem);
end

mpinv(mp) = 1:NumMod;                        % invert; model perm inverse

ir = 10000*mpinv(1:NumMod) + MNP;            % instance rank
                                             % nominal order within each model


[y,p] = sort(ir);                            % p is order for instances

q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];          % how to permute models, keep together
end



% -------------------------------------------- 

figure(31)
clf
if 0 > 1,
  T = RelSeq(p,q);                   % show forward and reversed models
  T(s+1,t+1) = 0;
  axis([1 2*s+1 1 s+1]);
else
  T = RelSeq(p,2*p-1);                  % show only forward models
  T(s+1,s+1) = 0;
  axis([1 s+1 1 s+1]);
end
pcolor(T);
hold on
axis ij
shading flat
caxis([-8 0]);
colormap('default')
map = colormap;
map = map(8:56,:);
%map = [1 1 1; map];
map = [map; 1 1 1];
colorbar('eastoutside');

ModelNames = Names(1:2:(t-1));

% ==========================================================================

figure(32)
clf

% ---------------------------------------- Graph individual instances

SS = S(p,q);                             % re-order scores acc. to perms

for m = 1:length(ModelNames),            % loop through models
  numthismodel = length(find(MN == MN(p(m))));  % number of instances in model
  [y,i] = sort(-SS(p(m),:));                    % scores
  plot(i(1:numthismodel),p(m)*ones(1,numthismodel),'.');
  hold on
end

axis([1 2*length(ModelNames) 1 length(ModelNames)]);
axis ij
title('For each sequence file, the top m scoring models, where m is the number of sequence files having the same model');

% ==========================================================================

figure(12)
clf

DD = abs(zDistance(max(-4,RelSeq)));   % how similarly sequences are scored
p = zOrderbySimilarity(DD);

q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];          % how to permute models, keep together
end

T = RelSeq(p,q);
T(s+1,t+1) = 0;
pcolor(T);
hold on
axis ij
shading flat
axis([1 2*s+1 1 s+1]);
caxis([-8 0]);
colormap('default')
map = colormap;
map = map(8:56,:);
%map = [1 1 1; map];
map = [map; 1 1 1];

for i = 1:s,
  plot([2*i 2*i+2 2*i+2 2*i 2*i]-1,[i i i+1 i+1 i],'k');
  hold on
end

colormap(map);
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',FullLab(p),'FontSize',5,'FontName','FixedWidth')
title('Scores of sequences against models and reversed models, relative to max score for those sequences');
xlabel('Models and reversed models, paired, in the order on the y axis');

figure(13)
clf

DD(end+1,end+1) = 0;
pcolor(DD(p,p));
hold on
axis ij
shading flat
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',FullLab(p),'FontSize',5,'FontName','FixedWidth')
title('Scores of sequences against models and reversed models, relative to max score for those sequences');
xlabel('Models and reversed models, paired, in the order on the y axis');

dt = strrep(datestr(now,31),':','_');

% ----------------------------------- order by how forward models scores seqs

figure(14)
clf

DD = abs(zDistance(max(-4,RelSeq(:,2:2:end)')));
p = zOrderbySimilarity(DD);

q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];          % how to permute models, keep together
end

T = RelSeq(p,q);
T(s+1,t+1) = 0;
pcolor(T);
hold on
axis ij
shading flat
axis([1 2*s+1 1 s+1]);
caxis([-8 0]);
colormap('default')
map = colormap;
map = map(8:56,:);
%map = [1 1 1; map];
map = [map; 1 1 1];

for i = 1:s,
  plot([2*i 2*i+2 2*i+2 2*i 2*i]-1,[i i i+1 i+1 i],'k');
  hold on
end

colormap(map);
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',FullLab(p),'FontSize',5,'FontName','FixedWidth')
title('Scores of sequences against models and reversed models, relative to max score for those sequences');
xlabel('Models and reversed models, paired, in the order on the y axis');

figure(10)
set(gcf,'renderer','painters')
set(gca,'fontsize',4)
orient landscape
saveas(gcf,['Sequence scores original order ' dt '_' num2str(counter) '_Correct.pdf'])

figure(12)
set(gcf,'renderer','painters')
set(gca,'fontsize',4)
orient landscape
saveas(gcf,['Sequence scores similarity order ' dt '_' num2str(counter) '_Correct.pdf'])

figure(13)
set(gcf,'renderer','painters')
set(gca,'fontsize',4)
orient landscape
saveas(gcf,['Distance between sequence scores ' dt '_' num2str(counter) '_Correct.pdf'])

