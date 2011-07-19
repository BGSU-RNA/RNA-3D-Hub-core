% pIsoScore4 returns a 4x4 matrix of probabilities depending on the specified
% interaction category Class and the specified Pair having this interaction.
% 85% of the probability goes to isosteric pairs (0 <= IDI <= 2), 10% to
% nearly isosteric pairs (2 < IDI <= 3.3), 4% to non, 1% to other

function [S] = pIsoScore4(Class,Code1,Code2,ExemplarIDI)

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

IDI = zExemplarIDI(Class,Code1,Code2,ExemplarIDI);

% ---------------------------------------- turn IDI into scores

for a = 1:4,
  for b = 1:4,
    if isnan(IDI(a,b)),
      IDI(a,b) = 100;                    % this pair is not possible
    end
  end
end

isosteric = sum(sum( IDI <= 2 ));        % how many pairs are isosteric
near      = sum(sum( (IDI > 2) .* (IDI <= 3.3)));  % how many near
non       = sum(sum( (IDI > 3.3) .* (IDI < 100))); % how many not
notallow  = 16 - isosteric - near - non;           % how many don't occur

% ---------------------------------------- normalize scores into probabilities

S = zeros(4,4);

if isosteric > 0,
  S = S + (0.85/isosteric) * (IDI <= 2);
end

if near > 0,
  S = S + (0.10/near) * (IDI > 2) .* (IDI <= 3.3);
end

if non > 0,
  S = S + (0.04/non) * (IDI > 3.3) .* (IDI < 100);
end

if notallow > 0,
  S = S + (0.01/notallow) * (IDI >= 100);
end

