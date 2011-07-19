
function [n] = xExtractNum(t)

c = 1;
a = '';
n = 0;

while c <= length(t),
  if any(t(c) == '0123456789'),
    a = [a t(c)];
    c = c + 1;
  elseif t(c) == ' ',
    c = length(t)+1;
  else
    c = c + 1;
  end
end

n = str2num(a);
