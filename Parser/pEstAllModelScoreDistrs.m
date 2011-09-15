% pEstAllModelScoreDistrs goes through all models of the current type, generates random sequence lengths, generates random sequences, scores them against the model itself, and saves the scores for use by other programs

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  break
end

top10 = 0;

sampSize = 100000;

ModelFNs = dir(['Models' filesep]);
keep = zeros(1,length(ModelFNs));
for m = 1:length(ModelFNs),
  if (length(ModelFNs(m).name) > 2),
    if (strcmp(ModelFNs(m).name(1:2),loopType)) && ~isempty(str2num(ModelFNs(m).name(4))),
      keep(m) = 1;
      ModelFNs(m).name = strrep(ModelFNs(m).name,'.txt','');
    end
  end
end
ModelFNs = ModelFNs(find(keep));

JAR3D_path;

%for m = 1:length(ModelFNs),
%for m = 1:120,
for m = 1:length(ModelFNs),
  MN = ModelFNs(m).name;
  BetterEmpDist(MN,loopType,sampSize,4,top10);
end

if top10 > 0,
  processor = 1;
  pCompareTopQuantiles
end
