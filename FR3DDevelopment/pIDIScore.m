% pIDIScore returns a 4x4 matrix of IDI's from the given basepair in the
% given family.

% Note to add:  AG water inserted is self-isosteric.  Call this class 14???

function [S] = pIDIScore(Edge,c,d)

if strcmp(class(Edge),'char'),
  Edge = xGetEdgeNums(Edge);
  Edge = Edge(1);
end

if strcmp(class(c),'char'),
  c = pLtoN(c);
end

if strcmp(class(d),'char'),
  d = pLtoN(d);
end

