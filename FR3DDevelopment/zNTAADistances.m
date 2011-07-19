
NTCenter = cat(1,File.NT.Center);

for aa = 1:length(File.AA),
  AACenter(aa,:) = File.AA(aa).Loc(2,:);    % use CA atom for each a.a.
end

DistanceNTAA = zDistance(NTCenter,AACenter);

DistanceNTAA = sparse( DistanceNTAA .* (DistanceNTAA < 20) );
DistanceNTAA = sparse( DistanceNTAA .* (DistanceNTAA < 10) );

spy(DistanceNTAA)

length(find(DistanceNTAA(127,:)))           % how many a.a.s are w/in 20 Ang

full(sum((DistanceNTAA > 0)'))

% min(min(zDistance(File.NT(127).Loc, File.AA(32).Loc)))
min(min(zDistance(File.NT(127).Fit, File.AA(32).Loc)))   % includes hydrogens
min(min(zDistance(File.NT(127).Sugar, File.AA(32).Loc))) % backbone

