% zLeadingZeros(n,m) converts integer n into a string and pads it with
% zeros to have total width m columns

function [T] = zLeadingZeros(n,m)

T = '';
for i = 1:m,
    T = [T '0'];
end

s = num2str(n);
if length(s) > m,
    T = s;
else
    T((end-length(s)+1):end) = s;
end
