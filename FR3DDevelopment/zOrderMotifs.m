% zOrderMotifs takes a set of motifs and matching matrix and orders the motifs by similarity to put the most similar motifs near each other in the list.  It takes as input a set of ModelNames; usually not all models are needed anyway.  It returns a permutation telling how to re-order ModelNames and a set of group-group distances in the same order as ModelNames

function [p,GDM,GSize] = zOrderMotifs(ModelNames,Sig,RSig,Method,Verbose)

if nargin < 2,
  Method = 1;
end

if nargin < 3,
  Verbose = 1;
end

loopType = ModelNames{1}(1:2);

if exist([loopType '_MatchingMatrix_Correct.mat']) == 2,
  load([loopType '_MatchingMatrix.mat']);               % load names, groups
  NumGroups = length(groups);

  for i = 1:NumGroups,
    GroupIndices{i} = find(ismember(names,groups{i}));
    GSize(i) = length(GroupIndices{i});

    if Verbose > 1,
      n1 = sort(groups{i});
      n2 = sort(names(GroupIndices{i}));

      for j = 1:length(GroupIndices{i}),
        fprintf('%-30s %-30s\n', n1{j}, n2{j});
      end
      fprintf('\n');
      pause
    end
  end

  MM_min = min(MM_init,MM_init');

  GroupDistI = zeros(NumGroups,NumGroups);
  GroupDistF = zeros(NumGroups,NumGroups);
  GroupDistM = zeros(NumGroups,NumGroups);

  for i = 1:NumGroups,
    for j = (i+1):NumGroups,
      GroupDistI(i,j) = median(reshape(MM_init(GroupIndices{i},GroupIndices{j}),1,[]));
      GroupDistI(j,i) = GroupDistI(i,j);

      GroupDistM(i,j) = median(reshape(MM_min(GroupIndices{i},GroupIndices{j}),1,[]));
      GroupDistM(j,i) = GroupDistM(i,j);

      GroupDistF(i,j) = median(reshape(MM_final(GroupIndices{i},GroupIndices{j}),1,[]));
      GroupDistF(j,i) = GroupDistF(i,j);
    end
  end

  for m = 1:length(ModelNames),

    GroupNumber(m) = str2num(ModelNames{m}(4:6));
  end

  GDI = GroupDistI(GroupNumber,GroupNumber);      % reduced-size distance matrix
  GDF = GroupDistF(GroupNumber,GroupNumber);      % reduced-size distance matrix
  GDM = GroupDistM(GroupNumber,GroupNumber);      % reduced-size distance matrix
  GSize = GSize(GroupNumber);

else
  Method = 3;
  GDM = zeros(length(ModelNames),length(ModelNames));
  GSize = 1:length(ModelNames);
end

% ---------------------------------------- Convert signatures to easier strings

for i = 1:length(Sig),
  m = Sig{i}(3:end);
%  m = strrep(m,'-','');
  m = strrep(m,'cWW','A');
  m = strrep(m,'tWW','B');
  m = strrep(m,'cWS','C');
  m = strrep(m,'cSW','D');
  m = strrep(m,'tWS','E');
  m = strrep(m,'tSW','F');
  m = strrep(m,'cWH','G');
  m = strrep(m,'cHW','H');
  m = strrep(m,'tWH','I');
  m = strrep(m,'tHW','J');
  m = strrep(m,'cHH','K');
  m = strrep(m,'tHH','L');
  m = strrep(m,'cHS','M');
  m = strrep(m,'cSH','N');
  m = strrep(m,'tHS','O');
  m = strrep(m,'tSH','P');
  m = strrep(m,'cSS','Q');
  m = strrep(m,'tSS','R');
  m = strrep(m,'bif','R');
  FSig{i} = m;

  m = RSig{i}(3:end);
%  m = strrep(m,'-','');
  m = strrep(m,'cWW','A');
  m = strrep(m,'tWW','B');
  m = strrep(m,'cWS','C');
  m = strrep(m,'cSW','D');
  m = strrep(m,'tWS','E');
  m = strrep(m,'tSW','F');
  m = strrep(m,'cWH','G');
  m = strrep(m,'cHW','H');
  m = strrep(m,'tWH','I');
  m = strrep(m,'tHW','J');
  m = strrep(m,'cHH','K');
  m = strrep(m,'tHH','L');
  m = strrep(m,'cHS','M');
  m = strrep(m,'cSH','N');
  m = strrep(m,'tHS','O');
  m = strrep(m,'tSH','P');
  m = strrep(m,'cSS','Q');
  m = strrep(m,'tSS','R');
  m = strrep(m,'bif','R');
  RSig{i} = m;
end

% ------------------------------------------ Calc distances between signatures
% ------------------------------------------ Align sequence signatures

FFDist = zeros(length(ModelNames));        % distances between forward sigs
FRDist = zeros(length(ModelNames));        % dist between forward and reversed

for i = 1:length(FSig),
  seq1 = FSig{i};                          % 
  for j = (i+1):length(FSig),
    seq2 = FSig{j};                        % 
    [matches,align1,align2,s1,s2] = dNeedlemanWunsch(seq1,seq2,0.99);
    FFDist(i,j) = abs(1 - (sum(s1==s2)) / length(s1));
    FFDist(j,i) = FFDist(i,j);

    seq2 = RSig{j};                        % 
    [matches,align1,align2,s1,s2] = dNeedlemanWunsch(seq1,seq2,0.99);
    FRDist(i,j) = abs(1 - (sum(s1==s2)) / length(s1));
    FRDist(j,i) = FRDist(i,j);
  end
end

SDist = min(FFDist,FRDist);                     % minimum distance

switch Method,
case 1,
  p = zOrderbySimilarity(GDM.*(0.5+SDist));    % mix signature and geometry
case 2,
  p = zOrderbySimilarity(GDM);                  % geometry only
case 3,
  p = zOrderbySimilarity(SDist);               % signature only
end

if Verbose > 0,
  ModelNames{p}
end

if Verbose > 1,
  figure(5)
  clf
  pcolor(GDI(p,p));
  hold on
  axis ij
  shading flat
  caxis([0 1]);
  colormap('default')
  map = colormap;
  map = map(56:-1:8,:);
  colormap(map);
  colorbar('eastoutside')
  set(gca,'YTick',(1:length(ModelNames))+0.5)
  set(gca,'YTickLabel',ModelNames(p),'FontSize',8,'FontName','FixedWidth')

  figure(6)
  clf
  pcolor(GDF(p,p));
  hold on
  axis ij
  shading flat
  caxis([0 1]);
  colormap('default')
  map = colormap;
  map = map(56:-1:8,:);
  colormap(map);
  colorbar('eastoutside')
  set(gca,'YTick',(1:length(ModelNames))+0.5)
  set(gca,'YTickLabel',ModelNames(p),'FontSize',8,'FontName','FixedWidth')

  figure(7)
  clf
  pcolor(GDM(p,p));
  hold on
  axis ij
  shading flat
  caxis([0 1]);
  colormap('default')
  map = colormap;
  map = map(56:-1:8,:);
  colormap(map);
  colorbar('eastoutside')
  set(gca,'YTick',(1:length(ModelNames))+0.5)
  set(gca,'YTickLabel',ModelNames(p),'FontSize',8,'FontName','FixedWidth')

  figure(8)
  clf
  pcolor(1-SDist(p,p));
  hold on
  axis ij
  shading flat
  caxis([0 1]);
  colormap('default')
  map = colormap;
  map = map(56:-1:8,:);
  colormap(map);
  colorbar('eastoutside')

  set(gca,'YTick',(1:length(ModelNames))+0.5)
  set(gca,'YTickLabel',ModelNames(p),'FontSize',8,'FontName','FixedWidth')

end

if Verbose > 1,
  [a,b] = find(GDM < 1);

  for i = 1:length(a),
    if a(i) > b(i),
      fprintf('%-60s and %-60s have distance %7.4f\n', ModelNames{a(i)}, ModelNames{b(i)}, GDF(a(i),b(i)));
    end
  end
end

