
function [n] = zSummarizeInteractions(File,Verbose)

for f = 1:length(File),
  E = fix(abs(File(f).Edge));
  C = (File(f).Crossing > 0);
  B = abs(File(f).BasePhosphate);
  n(f,1) = sum(sum( (1-C) .* (E == 1))) / 2;              % nested cWW
  n(f,2) = sum(sum( C .* (E == 1))) / 2;          % non-nested cWW
  n(f,3) = sum(sum( (1-C) .* (E > 1) .* (E < 13)))/2;     % nested non-cWW
  n(f,4) = sum(sum( C .* (E > 1) .* (E < 13)))/2; % non-nested non-cWW
  n(f,5) = sum(sum( (E > 19) .* (E < 25) ))/2;        % stacking
  n(f,6) = sum(sum( (B > 0) .* (B < 100)));           % base-phosphate
end

N{1} = 'nested cWW';
N{2} = 'non-nested cWW';
N{3} = 'nested non-cWW';
N{4} = 'non-nested non-cWW';
N{5} = 'stacking';
N{6} = 'base phosphate';

if Verbose > 0,
  fprintf('                      ');
  for f = 1:length(File),
    fprintf('%5s ', File(f).Filename);
  end
  fprintf('\n');
  for v = 1:6,
    fprintf('%20s ', N{v});
    for f = 1:length(File),
      fprintf('%5d ', n(f,v));
    end
    fprintf('\n');
  end
end
