
function [s] = pCoreSequence(a)

i = strfind(a,'*');
i = [i length(a)+1];

if length(i) > 1,
  s = a(2:(i(1)-2));              % up to first *
  for j = 1:(length(i)-1),
    s = [s '*' a((i(j)+2):(i(j+1)-2))];
  end
else
  s = a(2:(end-1));               % hairpin
end

