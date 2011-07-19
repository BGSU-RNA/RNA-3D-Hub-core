% zLogCounts does counts to make a stacked histogram colored by MinStructDist

function [logcounts,logbins] = zLogCounts(allscores,MinStructDist)

if nargin < 2,
  MinStructDist = ones(size(allscores));
end

  logscore = log(allscores);

  X = min(logscore):range(logscore)/30:max(logscore);

  clear logcounts
  clear logbins

  for d = 0:max(MinStructDist),
    i = find(MinStructDist == d);
    [logcounts(:,d+1),logbins] = hist(logscore(i),X);
  end
  
  [s,t] = size(logcounts);

  logcounts = [zeros(s,1) logcounts];

  logcounts = cumsum(logcounts')';

  logcounts = 2.2 * logcounts / max(max(logcounts));

  logcounts = 0.1*exp(logcounts);

%logcounts = logcounts/max(max(logcounts));
