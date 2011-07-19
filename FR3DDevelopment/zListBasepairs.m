
function [void] = zListBasepairs(File,start)

N = length(File.NT);

fprintf('[');

for i=start:N,
  for j=(i+1):N,
    a = fix(File.Edge(i,j));
    b = 0;
    switch a
      case    1,    b =  2;                     %  cWw
      case    2,    b =  3;                     %  tWw
      case    3,    b =  4;                     %  cWH
      case    4,    b =  5;                     %  tWH
      case    5,    b =  6;                     %  cWS
      case    6,    b =  7;                     %  tWS
      case    7,    b =  8;                     %  cHH
      case    8,    b =  9;                     %  tHH
      case    9,    b = 10;                     %  cHS
      case   10,    b = 11;                     %  tHS
      case   11,    b = 12;                     %  cSs
      case   12,    b = 13;                     %  tSs
      case   13,    b = 14;                     %  bif
      case   -1,    b =  2;                     %  cwW
      case   -2,    b =  3;                     %  twW
      case   -3,    b = 15;                     %  cHW
      case   -4,    b = 16;                     %  tHW
      case   -5,    b = 17;                     %  cSW
      case   -6,    b = 18;                     %  tSW
      case   -7,    b = 19;                     %  cHH
      case   -8,    b = 20;                     %  tHH
      case   -9,    b = 21;                     %  cSH
      case  -10,    b = 22;                     %  tSH
      case  -11,    b = 23;                     %  csS
      case  -12,    b = 24;                     %  tsS
      case  -13,    b = 25;                     %  bif
      case   21,    b =  0;                     %  s35
      case  -21,    b =  0;                     %  s53
      case   22,    b =  0;                     %  s33
      case  -22,    b =  0;                     %  s33
      case   23,    b =  0;                     %  s55
      case  -23,    b =  0;                     %  s55
      case  101,    b =  0;                     % ncWW
      case  102,    b =  0;                     % ntWW
      case  103,    b =  0;                     % ncWH
      case  104,    b =  0;                     % ntWH
      case  105,    b =  0;                     % ncWS
      case  106,    b =  0;                     % ntWS
      case  107,    b =  0;                     % ncHH
      case  108,    b =  0;                     % ntHH
      case  109,    b =  0;                     % ncHS
      case  110,    b =  0;                     % ntHS
      case  111,    b =  0;                     % ncSS
      case  112,    b =  0;                     % ntSS
      case -101,    b =  0;                     % ncWW
      case -102,    b =  0;                     % ntWW
      case -103,    b =  0;                     % ncHW
      case -104,    b =  0;                     % ntHW
      case -105,    b =  0;                     % ncSW
      case -106,    b =  0;                     % ntSW
      case -107,    b =  0;                     % ncHH
      case -108,    b =  0;                     % ntHH
      case -109,    b =  0;                     % ncSH
      case -110,    b =  0;                     % ntSH
      case -111,    b =  0;                     % ncSS
      case -112,    b =  0;                     % ntSS
      case -113,    b =  0;                     % nbif
      case  121,    b =  0;                     % ns35
      case -121,    b =  0;                     % ns53
      case  122,    b =  0;                     % ns33
      case -122,    b =  0;                     % ns33
      case  123,    b =  0;                     % ns55
      case -123,    b =  0;                     % ns55

    end
    if b > 0,
      fprintf('%d; ...\n', b);
%      fprintf('%s %s %s; ...\n', File.NT(i).Number, File.NT(j).Number, zBasepairCodeText(b));
    end
  end
end

fprintf('];\n');
