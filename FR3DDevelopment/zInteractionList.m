% zInteractionList(File,Display,Types) returns vectors indicating all the
% interactions in File
% Types is a 1 x 4 vector.  
% If Types(1) = 1, all basepairs will be included
% If Types(2) = 1, all stacking interactions will be included
% If Types(3) = 1, all near basepairs will be included
% If Types(4) = 1, all near stacking interactions will be included


function [Index1,Index2,InterCode] = zInteractionList(File,Display,Types,startindex,endindex)

N = length(File.NT);

if nargin < 5,
  endindex = N;
end

if nargin < 4,
  startindex = 1;
end

if nargin < 3,
  Types = [1 1 0 0];
end

if nargin < 2,
  Display = 0;
end

c = 1;                                          % interaction counter

for i=startindex:endindex,
  j = find(File.Edge(i,:));                     % columns of non-zero entries
  for k = 1:length(j),
    a = fix(File.Edge(i,j(k)));
    b = 0;
    if Types(1) > 0,
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
     end
    end

    if Types(2) > 0,
     switch a
      case   21,    b = 30;                     %  s35
      case  -21,    b = 31;                     %  s53
      case   22,    b = 32;                     %  s33
      case  -22,    b = 32;                     %  s33
      case   23,    b = 33;                     %  s55
      case  -23,    b = 33;                     %  s55
     end
    end

    if Types(3) > 0,
     switch a
      case  101,    b = 40;                     % ncWW
      case  102,    b = 41;                     % ntWW
      case  103,    b = 42;                     % ncWH
      case  104,    b = 43;                     % ntWH
      case  105,    b = 44;                     % ncWS
      case  106,    b = 45;                     % ntWS
      case  107,    b = 46;                     % ncHH
      case  108,    b = 47;                     % ntHH
      case  109,    b = 48;                     % ncHS
      case  110,    b = 49;                     % ntHS
      case  111,    b = 50;                     % ncSS
      case  112,    b = 51;                     % ntSS
      case  113,    b = 52;                     % nbif
      case -101,    b = 40;                     % ncWW
      case -102,    b = 41;                     % ntWW
      case -103,    b = 53;                     % ncHW
      case -104,    b = 54;                     % ntHW
      case -105,    b = 55;                     % ncSW
      case -106,    b = 56;                     % ntSW
      case -107,    b = 57;                     % ncHH
      case -108,    b = 58;                     % ntHH
      case -109,    b = 59;                     % ncSH
      case -110,    b = 60;                     % ntSH
      case -111,    b = 61;                     % ncSS
      case -112,    b = 62;                     % ntSS
      case -113,    b = 63;                     % nbif
     end
    end

    if Types(4) > 0,
     switch a
      case  121,    b = 70;                     % ns35
      case -121,    b = 71;                     % ns53
      case  122,    b = 72;                     % ns33
      case -122,    b = 72;                     % ns33
      case  123,    b = 73;                     % ns55
      case -123,    b = 73;                     % ns55
     end
    end

    if b > 0,
      Index1(c)    = i;
      Index2(c)    = j(k);
      InterCode(c) = b;
      c            = c + 1;
    end
  end
end

if Display > 0,
  for i = 1:length(Index1),
    fprintf('%5s %5s %4s\n', File.NT(Index1(i)).Number, File.NT(Index2(i)).Number, zInteractionCodeText(InterCode(i)));
  end
end