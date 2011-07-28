% zFitModifiedNucleotide(NT) identifies correspondence between atoms of modified nucleotide NT and its parent and uses those correspondences to find the rotation matrix that puts the modified base into standard orientation

function NT = zFitModifiedNucleotide(NT,Verbose)

if nargin < 2,
  Verbose = 0;
end

X = [[0 0 0]; [0 1 0]; [-1 0 0]];            % triangle in xy plane
Y = NT.Loc(1:3,:);                           % 3 points on modified base

% each entry of Corresp has two values, the atom name in the modified base, then the corresponding atom name in the parent base

Corresp = [];

[ModName,ModAtom,Parent,ParentAtom] = textread('ModifiedNucleotideAtomCorrespondence.txt','%s\t%s\t%s\t%s');

r = find(ismember(ModName,NT.Base));

if ~isempty(r),
  Parent = Parent{r(1)};
  Parent = find(Parent == 'ACGU');

  if ~isempty(Parent),

    MA = ModAtom(r);
    PA = ParentAtom(r);

    zStandardBases                % read in QM locations of atoms in 4 bases

    Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
    Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

    clear X Y

    for i = 1:length(MA),
      j = find(ismember(NT.AtomName,MA{i}));
      k = find(ismember(AtomNames(1:Lim(1,Parent),Parent),PA{i}));
      if isempty(j),
        fprintf('zFitModifiedNucleotide: %s is missing atom %s\n', NT.Base, MA{i});
      else
        X(i,:) = StandardLoc(k(1),:,Parent);  % ideal base atom locations
        Y(i,:) = NT.Loc(j(1),:);
      end
    end
  end
end

[Rot, sc, sh] = zBestTransformation(X,Y);      % find best rotation, shift

NT.Rot = Rot;                               % save the rotation matrix

NT.Fit = NT.Loc;                          % could do better if there is a QM-optimized version of this modified nucleotide; could also add hydrogens
NT.Syn = 0;

if Verbose > 0 && ~isempty(r),
  figure(83)
  clf
  zPlotStandardBase(Parent,1)
  clear VP
  VP.Sugar = 1;
  VP.AtOrigin = 1;
  clear F
  F.NT(1) = NT;
  zDisplayNT(F,1,VP);
  view(2)
end
