% zSubstituteNucleotides takes a given set of nucleotides and replaces individual nucleotides with different ones

% File.NT is a set of N nucleotides
% Codes is an Nx1 vector of codes, these are the new bases

function [Motif] = xSubstituteNucleotides(F,Indices,Codes,Verbose)

Motif = zSubFile(F,Indices);
N = length(Indices);

figure(1)
clf
VP.Sugar = 1;
zDisplayNT(Motif,'all',VP);

% define new base codes to use

if strcmp(class(Codes),'char'),
  Codes = 1*(Codes == 'A') + 2*(Codes == 'C') + 3*(Codes == 'G') + 4*(Codes == 'U');
end

% First, replace each Motif.NT.Fit with the new standard base

zStandardBases                % read in QM locations of atoms in 4 bases

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen
ACGU = 'ACGU';

for i = 1:N,                               % go through nucleotides

  dis = norm(Motif.NT(i).Fit(1,:) - Motif.NT(i).Sugar(1,:));

  if Motif.NT(i).Code ~= Codes(i),          % substitution is needed
    sh = Motif.NT(i).Fit(1,:)';            % current location of N1/N9 atom
    ce = Motif.NT(i).Center';               % current base center
    C  = Codes(i);                         % new code
    L2 = Lim(2,C);                         % num of atoms and hyd in new base
    X2 = StandardLoc(1:L2,:,C);            % ideal base & hyd locations
    sc = mean(X2);                         % center of standard base

    X2 = X2 - ones(L2,1)*sc;               % center standard base at center

%    X  = (sh*ones(1,L2) + Motif.NT(i).Rot*X2')';    % rotation of this base

                                           % align N1/N9 of old, new bases

    X  = (ce*ones(1,L2) + Motif.NT(i).Rot*X2')';    % rotation of this base

                                            % align center of old, new bases

    Motif.NT(i).Sugar = zRotateAboutPoint(Motif.NT(i).Sugar,Motif.NT(i).Sugar(13,:),sh,X(1,:),1)';

    Motif.NT(i).Fit(1:L2,:) = X;            % fitted locations of base, H
    Motif.NT(i).Code = C;
    Motif.NT(i).Base = ACGU(C);

    fprintf('Movement of N1/N9 by %7.4f Angstroms\n',norm(sh'-Motif.NT(i).Fit(1,:)));

  end

  fprintf('C1*-N1/N9 distance went from %7.4f to %7.4f\n', dis, norm(Motif.NT(i).Fit(1,:) - Motif.NT(i).Sugar(1,:)));

end

figure(2)
clf
zDisplayNT(Motif,'all',VP)

zLinkFigures(1:2)
