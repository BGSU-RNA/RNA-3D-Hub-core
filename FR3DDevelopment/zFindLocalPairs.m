
% File = zAddNTData('HighResolution_list',2);

for f = 1:length(File),
  LocalPair = zSimpleSecondaryStructure(File(f),10);
  save([pwd filesep 'PrecomputedData' filesep File(f).Filename '_LocalPair.mat'],'LocalPair');
end
