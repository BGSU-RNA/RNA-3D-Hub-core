% pMoreThanFlankingPairs(Names) goes through cell array Names, which must contain motif signatures, and identifies entries which consist of more than just flanking cWW pairs

function [MoreThan] = pMoreThanFlankingPairs(Names,LoopType)

if nargin < 2,
  LoopType = 'IL';
end

switch LoopType
case 'HL'
  NumcWW = 1;
case 'IL'
  NumcWW = 2;
otherwise
  NumcWW = 3;
end

MoreThan = [];
for m = 1:length(Names),
  MoreThan(m) = 0;
  if ~isempty(strfind(Names{m},'-t')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},'-cS')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},'-cH')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},'-cWS')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},'-cWH')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},',t')),
    MoreThan(m) = 1;
  elseif ~isempty(strfind(Names{m},',c')),
    MoreThan(m) = 1;
  elseif length(strfind(Names{m},'cWW')) > 2,
    MoreThan(m) = 1;
  end
end
