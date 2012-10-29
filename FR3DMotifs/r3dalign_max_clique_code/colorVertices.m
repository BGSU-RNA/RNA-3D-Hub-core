%Greedy coloring algorithm; G is the edge matrix where zeros denote
%adjacent vertices

function [Colors] = colorVertices(G)

Colors = zeros(1,length(G));
numColors=0;
S=sum(G);
F=find(S==0);
GG=G;
U=setdiff(1:length(G),F);
G=GG(U,U);
tColors=zeros(1,length(G));
for i = 1:length(G)     
%      N=find(G(i,1:i-1)==0);     %get the already colored neighbors
%      J=sort(tColors(N));
   J=sort(tColors(G(i,1:i-1)==0));
   D=diff(J);
   L=find(D>1);
   if ~isempty(L)
      NC=J(L(1))+1;
      tColors(i) = NC;
   else
      if ~isempty(J)
         NC=J(end)+1;
      else
         NC=1;
      end
      if NC > numColors
         numColors=numColors+1;
         tColors(i)=numColors;
      else
         tColors(i) = NC;
      end
   end
end
Colors(U)=tColors+length(F);
Colors(F)=1:length(F);
