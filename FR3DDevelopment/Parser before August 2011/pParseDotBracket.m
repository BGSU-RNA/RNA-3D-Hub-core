% pParseDotBracket(DB,LeftColumn) identifies the columns in which each basepair occurs.  It assumes that the parentheses are nested.  LeftColumn tells where in the original sequence the first (left-most) element of the current string belongs

function [PairColumn] = pParseDotBracket(DB,LeftColumn)

fprintf('%4d %s\n', LeftColumn, DB);

i = min(find(DB == '('));                % location of first (

if isempty(i),                           % no parenthesis here
  PairColumn = [];                       % return an empty list of pairs
else
  p = 1*(DB(i:end) == '(') - 1*(DB(i:end) == ')');
                                         % 1 and -1 for ( and ), respectively
  q = cumsum(p);
  r = i-1+min(find(q == 0));             % first time parentheses balance out

  if isempty(r),
    fprintf('pParseDotBracket did not find a matching parenthesis\n');
    fprintf('%s\n', DB);
  else
    PairColumn = LeftColumn - 1 + [i r];   % matched pair of columns
    PairColumn = [PairColumn; pParseDotBracket(DB(i+1:r-1),LeftColumn+i)];
    PairColumn = [PairColumn; pParseDotBracket(DB((r+1):end),LeftColumn+r)];  
  end
end
