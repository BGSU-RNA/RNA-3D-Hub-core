% zRotateAboutPoint(x,p,a,b,s) rotates points in the 3xN matrix x about point p in such a way that point a goes to the line from p to b.  If s > 0, the points x are then scaled so that a maps to b.  The rotated points are returned in y.

function [y] = zRotateAboutPoint(x,p,a,b,s)

[M,N] = size(x);
if M ~= 3,
  x = x';                                 % hope for the best
  [M,N] = size(x);
end

if nargin < 5,
  s = 0;
end

p = reshape(p,3,1);                       % make sure p is 3x1
a = reshape(a,3,1);                       % make sure a is 3x1
b = reshape(b,3,1);                       % make sure b is 3x1

aa = (a-p) / norm(a-p);
bb = (b-p) / norm(b-p);

z = cross(aa,bb);                         % axis about which to rotate
z = z / norm(z);
n = cross(z,bb);                          % perpendicular to the axis
n = n / norm(n);

t = acos(aa'*bb);                         % angle to rotate through

if s > 0,
  f = norm(b-p)/norm(a-p);                % factor to stretch/shrink
else
  f = 1;                                  % no stretch/shrink
end

R = [bb n z];                             % rotate from e_1 e_2 e_3 to bb n z
A = [ f*cos(t) -f*sin(t) 0; f*sin(t) f*cos(t) 0; 0 0 1];  % rotate by angle t

y = p*ones(1,N) + R * A * R' * (x-p*ones(1,N));

if 0 > 1,
aa
bb
z
n

R
A
end


if 0 > 1,
  zRotateAboutPoint([1 0 1; 3 0 0; 1 0 0]',[0 0 0],[1 0 0],[0 1 0])
  zRotateAboutPoint([1 0 1; 3 0 0; 1 0 0]',[0 0 0],[1 0 0],[0 2 0],1)
  zRotateAboutPoint([1 0 1; 3 0 0; 1 0 0]',[0 0 0],[1 0 0],[0 0.5 2],1)
  zRotateAboutPoint([1 0 0; 3 0 0; 0 1 0; 0 0 1]',[0 0 0],[1 0 0],[1 1 1],0)
  zRotateAboutPoint([1 0 0; 3 0 0; 0 1 0; 0 0 1]',[1 1 0],[1 0 0],[0 1 0],0)

  zRotateAboutPoint([1 0 0; 3 0 0; 0 1 0; 0 0 1]',[0 0 0],[1 0 0],[0 0 1],0)
end
