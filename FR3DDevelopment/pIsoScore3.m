% pIsoScore3 returns a 4x4 matrix of scores depending on the specified
% interaction category Class and the specified Pair having this interaction.

function [S] = pIsoScore3(Class,Code1,Code2,ExemplarIDI,ExemplarFreq)

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

S = 1 ./ (1 + (0.5*IDI).^2);                   % score decreases with IDI

S = max(S,0.05);                         % set a minimum score for all pairs

% ---------------------------------------- multiply scores by frequencies

F = ExemplarFreq{fix(abs(Class))};

F = F / sum(sum(F));                     % normalize frequencies
F = max(F,0.01);                         % make minimum frequency about 1%

S = S .* F;

% ---------------------------------------- normalize scores into probabilities

S = S / (sum(sum(S)));                   % normalize
