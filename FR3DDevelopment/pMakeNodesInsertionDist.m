% pMakeNodesInsertionDist(codes.Type) creates a length and letter distribution from the given codes, according to the Type of location where the insertions occur

function [ID] = pMakeNodesInsertionDist(codes,Type,ZeroProb)

if nargin < 2,
  Type = 'Unknown';
end

L = length(codes) + 1;              % shift by one; entry 1 is for 0 insertions

D = [0.9 0.1];                          % default

switch Type
case {'Basepair','Unknown'}
  a = 1:(L+2);
  D = (0.1).^(abs(a-L));                % set the length distribution
case 'Cluster'
  a = 1:(L+1);
  D = (0.1).^(abs(a-L));                % set the length distribution
end

ID.LengthDist = D/sum(D);               % normalize the length distribution

if length(codes) == 0 && strcmp(Type,'Basepair'),
  ID.LengthDist = [0.995 0.0045 0.0005]; % probability of a new insertion
end

if length(codes) == 0 && strcmp(Type,'Cluster'),
  ID.LengthDist = [0.999 0.0009 0.0001]; % probability of a new insertion
end

switch Type
case {'Basepair','Unknown'}
  LettPref = 1;                         % strength of preference for observed
case 'Cluster'
  LettPref = 2;                         % stronger preference
case 'Initial'
  LettPref = 0.5;
end

LD = [1 1 1 1]';                      % default inserted letter distribution

for c = 1:length(codes),
  LD(codes(c)) = LD(codes(c))+LettPref; % give preference to observed letters
end  

ID.LetterDist = LD/sum(LD);             % normalize the letter preference

if 1 > 1,                               % display results
  Lett = 'ACGU';
  codes
  Lett(codes)
  ID.LengthDist
  ID.LetterDist
end
