
function [j] = zInvertFunction(i)

n = length(i);
m = max(i);

for a = n:-1:1,
  if i(a) > 0,
    j(i(a)) = a;
  end
end

z = find(j == 0);
j(z) = (n+1) * ones(size(z));