
function [F] = pModifyEdge(F,N1,N2,EdgeText)

if strcmp(class(EdgeText),'char'),
  e = zEdgeFromEdgeText(EdgeText);
else
  e = EdgeText;
end

i1 = zIndexLookup(F,{N1});
i2 = zIndexLookup(F,{N2});

F.Edge(i1,i2) =  e;
F.Edge(i2,i1) = -e;

