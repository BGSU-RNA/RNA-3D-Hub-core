
function [Distance, Score] = zBasePhosphateGeometry(BPh)

if 1 <= BPh && BPh <= 19,

zStandardBases

for Code = 1:4,
    switch Code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

  for w = 1:length(h),
    HLoc(e(w),:) = StandardLoc(h(w),:,Code);
  end
end

BPhCodes{1} = [1 2 3 4];                     % codes for A
BPhCodes{2} = [5 6 7 8 9 18];                % codes for C
BPhCodes{3} = [10 11 12 13 14 19];           % codes for G
BPhCodes{4} = [15 16 17];                    % codes for U

HLoc( 7,:) = (HLoc( 6,:)+HLoc( 8,:))/2;      % average these
HLoc(18,:) = (HLoc( 6,:)+HLoc( 8,:))/2;      % average these
HLoc(12,:) = (HLoc(11,:)+HLoc(13,:))/2;      % average these
HLoc(19,:) = (HLoc(11,:)+HLoc(13,:))/2;      % average these

BPhDist = zDistance(HLoc,HLoc);

for c = 1:4,
  Distance(c) = min(BPhDist(BPh,BPhCodes{c}));
end

Score = 1./(1+2*Distance.^2);                  % convert to probabilities
Score = Score / sum(Score);                    % normalize

Lett = 'ACGU';

if 0 > 1,
  BPhCodes{1} = [1 2 3 4];                     % codes for A
  BPhCodes{2} = [5 6 7 8 9 18];                % codes for C
  BPhCodes{3} = [10 11 12 13 14 19];           % codes for G
  BPhCodes{4} = [15 16 17];                    % codes for U
  for BPh = 1:19,
    for c = 1:4,
      if any(BPh == BPhCodes{c}),
        Letter = Lett(c);
      end
    end
    [D,S] = zBasePhosphateGeometry(BPh);
    fprintf('BPh code %d is %s made by %s\n', BPh, zBasePhosphateText(BPh), Letter);
    D
    S
  end
end

else

Distance = [];
Score = [];

end
