% G is the edge matrix of a graph
% adjacent vertices are represented by zeros in G
function [Clique1] = rGetMaximalClique(G)

Clique1 = [];
S=1:length(G);
while isempty(S) == false 
   T=sum(G);
   if sum(sum(T))==0    
      Clique1 = horzcat(Clique1,S);
      break;
   else
      I=find(T==0);
      if ~isempty(I)
         v=S(I);
         Clique1 = horzcat(Clique1,v);
         Z=setdiff(1:length(S),I);
         G=G(Z,Z);
         S(I)=[];
      else  
         [C I] = min(T);                    % C is value of smallest sum of column and I is the column index into S
         v=S(I);                            % vertex to add to clique
         Clique1 = horzcat(Clique1,v);      % add that vertex to clique 
         Z=[I find(G(I,:))];                % Z contains the newly added vertex and its non-adjacent vertices
         Q=setdiff(1:length(S),Z);          % Remove Z vertices from list of eligible vertices for clique
         G=G(Q,Q);
         S(Z)=[];
      end
   end
end 