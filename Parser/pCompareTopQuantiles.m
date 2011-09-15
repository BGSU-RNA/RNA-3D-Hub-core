% pCompareTopQuantiles parses the top 10th percentile sequences from each model against all models and compares them

if ~exist('loopType'),
  loopType = 'HL';
  loopType = 'IL';
end

if 0>1,
  loopType = 'IL';
  processor = 1;
  load(['pCompareTopQuantiles' num2str(processor)]);
end

% --------------------------------- load model names

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

% ---------------------------------- Remove models with no basepair beside
% ---------------------------------- flanking cWW

if basepaironly > 0,
  Keep = pMoreThanFlankingPairs(ModelNames,loopType);

  i = find(Keep);                  % indices of models to keep

  fprintf('Keeping %d models out of %d because they have more than just flanking cWW pairs\n', length(i), length(ModelNames));

  ModelNames = ModelNames(i);
end

% --------------------------------- sort by increasing model size

clear msize

for i = 1:length(ModelNames),
  tem = ModelNames{i};
  if isempty(strfind(tem,'Helix')),
    msize(i) = str2num(tem(8:9));
  else
    msize(i) = str2num(tem(14:15));
  end
end

[y,i] = sort(msize);
ModelNames = ModelNames(i);                % sort by decreasing model size

% --------------------------------- parse these top sequences against all

%          1  2  3  4  5  6  7  8  9 10 11   12 13   14   15   16    17     18
cutoffs = [0 90 91 92 93 94 95 96 97 98 98.5 99 99.5 99.8 99.9 99.99 99.999 100]/100;

Counts{length(ModelNames),length(ModelNames)} = [];

JAR3D_path

L = [1 120; 121 240; 241 356];

didone = 0;

for i = processor:3:length(ModelNames),           % loop through all models
  SFN = ['Top10_' ModelNames{i}(1:6) '.fasta'];

  MN = ModelNames{i};
  Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,MN);
  MyQuant = webJAR3D.getQuantilesB(Scores, MN(4:6), loopType);

  for j = 1:length(ModelNames),
   onemodel = tic;
   if isempty(Counts{i,j}),
    if j ~= i,

      didone = 1;

      NewMN = ModelNames{j};
      fprintf('(%d,%d) Sequences of %s against %s\n', i, j, MN, NewMN);

      Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,NewMN);
      quan = tic;
      Quant = webJAR3D.getQuantilesB(Scores, NewMN(4:6), loopType);
      toc(quan)

if 0>1,
      clf
      plot(MyQuant,Quant,'.');
      drawnow
      pause
end
    else
      NewMN = MN;
      Quant = MyQuant;
    end

    for a = 1:length(cutoffs),
      for b = 1:length(cutoffs),
        Counts{i,j}(a,b) = length(find((MyQuant >= cutoffs(a)) .* (Quant >= cutoffs(b))));
      end
    end

    if Counts{i,j}(2,2) == 0,           % no intersection
      Counts{j,i} = zeros(length(cutoffs));      % don't bother checking here
    end

    Counts{i,j}([1 2 7 10 12 15], [1 2 7 10 12 15])

    toc(onemodel)
   end
  end

  if didone > 0,
    save(['pCompareTopQuantiles' num2str(processor)],'Counts','ModelNames');
    whos

    w = 0;
    for i = 1:length(ModelNames),
      for j = 1:length(ModelNames),
        if ~isempty(Counts{i,j}),
          w = w + 1;
        end
      end
    end
    fprintf('%d of %d entries done before updating\n',w, length(ModelNames)^2);

    CC = Counts;                                   % save counts

    tim = clock;
    while mod(fix(tim(6)),10) ~= 3*processor,       % separate times
      tim = clock;
    end
  
    for p = 1:3,
      if p ~= processor && exist(['pCompareTopQuantiles' num2str(p) '.mat']) == 2,
        load(['pCompareTopQuantiles' num2str(p)]);    % careful!  clears Counts
        for i = 1:length(ModelNames),
          for j = 1:length(ModelNames),
            if isempty(CC{i,j}) && ~isempty(Counts{i,j}),
              CC(i,j) = Counts(i,j);
            end
          end
        end
      end
    end

    Counts = CC;

    w = 0;
    for i = 1:length(ModelNames),
      for j = 1:length(ModelNames),
        if ~isempty(Counts{i,j}),
          w = w + 1;
        end
      end
    end
    fprintf('%d of %d entries done after updating\n',w, length(ModelNames)^2);
  end
end

% -------------------------------------------------------------------------

if 0 > 1,

  CC = Counts;

  C = Counts;

%          1  2  3  4  5  6  7  8  9 10 11   12 13   14   15   16    17     18
cutoffs = [0 90 91 92 93 94 95 96 97 98 98.5 99 99.5 99.8 99.9 99.99 99.999 100]/100;

  if ~exist('ModelNames'),
    for i = 1:length(ModelFNs),
      ModelNames{i} = ModelFNs(i).name;
    end
  end

  s = length(ModelNames);

  Overlap = zeros(s,s);
  PercOver = zeros(s,s);
  MaxOverPercent = zeros(s,s);

  cu = 2;                                              % cutoff to use

  % Overlap(i,j) is the percentage overlap of cutoffs(cu)
  % PercOver(i,j) is the highest percentile above which there is overlap
  % MaxOverPercent(i,j) is the largest percentage overlap above 90%, 95%, etc.
  for i = 1:s,
    for j = 1:s,
      if isempty(C{i,j}),
        Overlap(i,j) = NaN;
        PercOver(i,j) = NaN;
        MaxOverPercent(i,j) = NaN;
      elseif length(C{i,j}(:,1)) == 1,
        Overlap(i,j)  = 0;
        PercOver(i,j) = 0;
        MaxOverPercent(i,j) = 0;
      else
        MaxOverPercent(i,j) = 0;
        Overlap(i,j) = C{i,j}(cu,cu) / C{i,j}(cu,1);
wroteone = 0;
        for c = 2:17,
          if C{i,j}(c,c) > 0,
            PercOver(i,j) = cutoffs(c);
fprintf('%7.4f ', C{i,j}(c,c) / C{i,j}(c,1));
wroteone = 1;
            MaxOverPercent(i,j) = max(MaxOverPercent(i,j),C{i,j}(c,c)/C{i,j}(c,1));
          end
        end
if wroteone > 0,
  fprintf('\n');
end
      end
    end
  end

  for i = 1:s,
    MNLab{i} = ModelNames{i}(4:(end-4));
    MNLab{i} = strrep(MNLab{i},'----','-');
    MNLab{i} = strrep(MNLab{i},'---','-');
    MNLab{i} = strrep(MNLab{i},'--','-');
    MNLab{i} = strrep(MNLab{i},'-','-');
  end

per = 1:length(ModelNames);

  figure(6)
  clf
  MinOver = Overlap;
  MinOver = min(Overlap,Overlap');
  per = zOrderbySimilarity(max(max(MinOver))-MinOver);
  MinOver(end+1,end+1) = 0;
  pcolor(MinOver(per,per));
  hold on
  axis ij
  shading flat
%  caxis([0 0.1]);
  colormap('default')
  map = colormap;
%  map = map(56:-1:8,:);
  map = map(8:56,:);
  colormap(map);
  colorbar('eastoutside')
%  set(gca,'XTick',(1:length(MNLab))+0.5)
%  set(gca,'XTickLabel',MNLab(per),'FontSize',8,'FontName','FixedWidth')
  set(gca,'YTick',(1:6:length(MNLab))+0.5)
  set(gca,'YTickLabel',MNLab(per(1:6:length(MNLab))),'FontSize',6,'FontName','FixedWidth')
  title(['Overlap of ' num2str(100*cutoffs(cu)) 'th percentiles']);

  figure(7)
  clf
  PercOver = min(PercOver,PercOver');           
  per = zOrderbySimilarity(max(max(PercOver))-PercOver);
  PercOver(end+1,end+1) = 0;
  pcolor(100*PercOver(per,per));
  hold on
  axis ij
  shading flat
  caxis([85 100]);
  colormap('default')
  map = colormap;
%  map = map(56:-1:8,:);
  map = map(8:56,:);
  colormap(map);
  colorbar('eastoutside')
%  set(gca,'XTick',(1:length(MNLab))+0.5)
%  set(gca,'XTickLabel',MNLab(per),'FontSize',8,'FontName','FixedWidth')
  set(gca,'YTick',(1:6:length(MNLab))+0.5)
  set(gca,'YTickLabel',MNLab(per(1:6:length(MNLab))),'FontSize',6,'FontName','FixedWidth')
  title(['Non-empty overlap of percentiles']);

  figure(8)
  clf
  MinMax = min(MaxOverPercent,MaxOverPercent');
  per = zOrderbySimilarity(max(max(MinMax))-MinMax);
  MinMax(end+1,end+1) = 0;
  pcolor(100*MinMax(per,per));
  hold on
  axis ij
  shading flat
  colormap('default')
  map = colormap;
  map = map(8:56,:);
  colormap(map);
  colorbar('eastoutside')
%  set(gca,'XTick',(1:length(MNLab))+0.5)
%  set(gca,'XTickLabel',MNLab(per),'FontSize',8,'FontName','FixedWidth')
  set(gca,'YTick',(1:6:length(MNLab))+0.5)
  set(gca,'YTickLabel',MNLab(per(1:6:length(MNLab))),'FontSize',6,'FontName','FixedWidth')
%  set(gca,'YTick',(1:length(MNLab))+0.5)
%  set(gca,'YTickLabel',MNLab(per(1:length(MNLab))),'FontSize',6,'FontName','FixedWidth')
  title(['Highest percentage overlap of top-ranking sequences'],'FontSize',14,'FontName','VariableWidth');

if 10 > 1,
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

  for i = 1:length(ModelNames),
    j = per(i);

    fprintf('%4d\t%sGroup_%s.html\t%s\n', i, BaseURL, ModelNames{j}(4:6), ModelNames{j});


    if i < length(ModelNames),
%    MinMax(per(i),per(i+1))
      if MinMax(per(i),per(i+1)) < 0.20,
        fprintf('\n');
      end
    end
  end
end

  Counts = CC;

end
