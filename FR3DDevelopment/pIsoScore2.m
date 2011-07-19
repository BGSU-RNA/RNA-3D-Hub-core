% pIsoScore2 returns a 4x4 matrix of scores depending on the specified
% interaction category Class and the specified Pair having this interaction.

function [S,IDI] = pIsoScore2(Class,Code1,Code2,ExemplarIDI)

if nargin < 4,
  load PairExemplars
end

if strcmp(class(Class),'char'),
  Class = xGetEdgeNums(Class);
end

if strcmp(class(Code1),'char'),
  Code1 = pLtoN(Code1);
end

if strcmp(class(Code2),'char'),
  Code2 = pLtoN(Code2);
end

% ---------------------------------------- look up IDI values

% [Class Code1 Code2]

IDI = zExemplarIDI(Class,Code1,Code2,ExemplarIDI);

% ---------------------------------------- turn IDI into scores

for a = 1:4,
  for b = 1:4,
    if isnan(IDI(a,b)),
      IDI(a,b) = 100;                    % this pair is not possible
    end
  end
end

%S = 1 * (IDI < 2.0) + 1 * exp(-0.2*(IDI-2.0)) .* (IDI >= 2.0);

S = 1 ./ (1 + (0.5*IDI).^2);                   % score decreases with IDI

S = max(S,0.05);                         % set a minimum score for all pairs

% ---------------------------------------- normalize scores into probabilities

S = S / (sum(sum(S)));                   % normalize
