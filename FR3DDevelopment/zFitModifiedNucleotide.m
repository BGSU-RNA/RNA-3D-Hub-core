% zFitModifiedNucleotide(NT) identifies correspondence between atoms of modified nucleotide NT and its parent and uses those correspondences to find the rotation matrix that puts the modified base into standard orientation

function NT = zFitModifiedNucleotide(NT,Verbose)

if nargin < 2,
  Verbose = 0;
end

X = [[0 0 0]; [0 1 0]; [-1 0 0]];            % triangle in xy plane
Y = NT.Loc(1:3,:);                           % 3 points on modified base

% each entry of Corresp has two values, the atom name in the modified base, then the corresponding atom name in the parent base

Corresp = [];

switch NT.Base,

case '1MA',
  Parent = 1;                   % adenosine is the parent
  Corresp{1}  = {'N9','N9'};
  Corresp{2}  = {'C4','C4'};
  Corresp{3}  = {'C5','C5'};
  Corresp{4}  = {'N7','N7'};
  Corresp{5}  = {'C8','C8'};
  Corresp{6}  = {'N3','N3'};
  Corresp{7}  = {'C2','C2'};
  Corresp{8}  = {'N1','N1'};
  Corresp{9}  = {'C6','C6'};
  Corresp{10} = {'N6','N6'};

end


if ~isempty(Corresp),
  zStandardBases                % read in QM locations of atoms in 4 bases

  Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  clear X Y

  for i = 1:length(Corresp),
    j = find(ismember(NT.AtomName,Corresp{i}{1}));
    k = find(ismember(AtomNames(1:Lim(1,Parent),Parent),Corresp{i}{2}));
    if isempty(j),
      fprintf('zFitModifiedNucleotide: %s is missing atom %s\n', NT.Base, Corresp{i}{1});
    else
      X(i,:) = StandardLoc(k(1),:,Parent);  % ideal base atom locations
      Y(i,:) = NT.Loc(j(1),:);
    end
  end

end

[r, sc, sh] = zBestTransformation(X,Y);      % find best rotation, shift

NT.Rot = r;                               % save the rotation matrix

NT.Fit = NT.Loc;                          % could do better if there is a QM-optimized version of this modified nucleotide; could also add hydrogens
NT.Syn = 0;

if Verbose > 0 && ~isempty(Corresp),
  figure(83)
  clf
  zPlotStandardBase(Parent,1)
  VP.Sugar = 1;
  VP.AtOrigin = 1;
  F.NT(1) = NT;
  zDisplayNT(F,1,VP);
  view(2)
end
