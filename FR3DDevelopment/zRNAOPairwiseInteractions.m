% zRNAOPairwiseInteractions(e) converts internal codes for pair interactions into text

% e is the internal code

function [E] = zRNAOPairwiseInteractions(e,C1,C2)

if nargin < 1,
  e = 1;
end

E = [];

for i=1:length(e),
  switch fix(e(i))
    case    1,    E = [E 'pairs_with_CWW'];
    case    2,    E = [E 'pairs_with_TWW'];
    case    3,    E = [E 'pairs_with_CWH'];
    case    4,    E = [E 'pairs_with_TWH'];
    case    5,    E = [E 'pairs_with_CWS'];
    case    6,    E = [E 'pairs_with_TWS'];
    case    7,    E = [E 'pairs_with_CHH'];
    case    8,    E = [E 'pairs_with_THH'];
    case    9,    E = [E 'pairs_with_CHS'];
    case   10,    E = [E 'pairs_with_THS'];
    case   11,    E = [E 'pairs_with_CSS'];
    case   12,    E = [E 'pairs_with_TSS'];
    case   -1,    E = [E 'pairs_with_CWW'];
    case   -2,    E = [E 'pairs_with_TWW'];
    case   -3,    E = [E 'pairs_with_CHW'];
    case   -4,    E = [E 'pairs_with_THW'];
    case   -5,    E = [E 'pairs_with_CSW'];
    case   -6,    E = [E 'pairs_with_TSW'];
    case   -7,    E = [E 'pairs_with_CHH'];
    case   -8,    E = [E 'pairs_with_THH'];
    case   -9,    E = [E 'pairs_with_CSH'];
    case  -10,    E = [E 'pairs_with_TSH'];
    case  -11,    E = [E 'pairs_with_CSS'];
    case  -12,    E = [E 'pairs_with_TSS'];
    case   21,    E = [E 'stacks_three_prime_face_five_prime_face'];
    case  -21,    E = [E 'stacks_five_prime_face_three_prime_face'];
    case   22,    E = [E 'stacks_three_prime_face_three_prime_face'];
    case  -22,    E = [E 'stacks_three_prime_face_three_prime_face'];
    case   23,    E = [E 'stacks_five_prime_face_five_prime_face'];
    case  -23,    E = [E 'stacks_five_prime_face_five_prime_face'];
    case  101,    E = [E 'ncWw'];
    case  102,    E = [E 'ntWw'];
    case  103,    E = [E 'ncWH'];
    case  104,    E = [E 'ntWH'];
    case  105,    E = [E 'ncWS'];
    case  106,    E = [E 'ntWS'];
    case  107,    E = [E 'ncHh'];
    case  108,    E = [E 'ntHh'];
    case  109,    E = [E 'ncHS'];
    case  110,    E = [E 'ntHS'];
    case  111,    E = [E 'ncSs'];
    case  112,    E = [E 'ntSs'];
    case  113,    E = [E 'nbif'];
    case  114,    E = [E 'nRib'];
    case -101,    E = [E 'ncwW'];
    case -102,    E = [E 'ntwW'];
    case -103,    E = [E 'ncHW'];
    case -104,    E = [E 'ntHW'];
    case -105,    E = [E 'ncSW'];
    case -106,    E = [E 'ntSW'];
    case -107,    E = [E 'nchH'];
    case -108,    E = [E 'nthH'];
    case -109,    E = [E 'ncSH'];
    case -110,    E = [E 'ntSH'];
    case -111,    E = [E 'ncsS'];
    case -112,    E = [E 'ntsS'];
    case -113,    E = [E 'nbif'];
    case -114,    E = [E 'nRib'];
    case  121,    E = [E 'ns35'];
    case -121,    E = [E 'ns53'];
    case  122,    E = [E 'ns33'];
    case -122,    E = [E 'ns33'];
    case  123,    E = [E 'ns55'];
    case -123,    E = [E 'ns55'];
    otherwise     E = [E '----'];
  end


return

EE = E;
EE((end-1):end) = upper(E((end-1):end));

if nargin == 4,
  if any(fix(abs(e)) == [11 12 111 112]),
    EE = E;                                     % reset to original
  elseif any(abs(e) == [1 7 8]) && (C1 == C2),
    EE = E;
  end
end

if nargin == 3,
  if any(fix(abs(e)) == [11 12 111 112]),
    EE = E;                                     % reset to original
  elseif any(abs(e) == [1 7 8]) && (any(C1 == [1 6 11 16])),
    EE = E;
  end
end

if any(fix(abs(e)) == [13 14 113 114]),
  EE = E;
end

E = EE;

if abs(e(i)) < 29,                     % pairing or stacking
  if (fix(e(i)) == e(i)) || (Detail == 0),
    E = [E ' '];
  else
    d = abs(e(i)) - fix(abs(e(i)));
    d = fix(d*10+0.0001);
    d = max(1,d);
    T = 'abcdefghijkl';
    E = [E T(d)];
  end
 end
end

