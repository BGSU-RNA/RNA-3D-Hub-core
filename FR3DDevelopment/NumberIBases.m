%----------------------------------------------------------------
function Node=NumberIBases(Node)

L = length(Node.IBases(:,1));
[a,b] = unique(Node.IBases);
for k=1:length(a)
   [c,d]= find(Node.IBases==a(k));
   for t=1:length(c)
        IBases(c(t),d(t))=k;
        IBasesOld(c(t),d(t))=a(k);
   end
end
Node.IBases=IBases;
Node.IBasesOld=IBasesOld;
