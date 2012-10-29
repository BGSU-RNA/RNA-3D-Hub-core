function [CBC] = rCliqueBB(G)

R=rReorderVerts(G,'ascend');
G=G(R,R);
V = 1:length(G);  % V = [ 1 2 3 ... n];
L=[];         %L will be a cell array that serves as a stack
              %L{i,1} will contain Vi
              %L{i,2} will contain Ui
              %L{i,3} will contain coloring of obtained by using colors
              %of branch above;  its length is Xhat
StackCt=1;                
S=sum(G);
F=find(S==0);
L{StackCt,1}=setdiff(V,F);
L{StackCt,2}=F;

CBC=rGetMaximalClique(G);
CBCsize=length(CBC);
 
C=colorVertices(G);
L{StackCt,3}=C;

while StackCt > 0
     
   Vt=L{StackCt,1};
   Ut=L{StackCt,2};
     
   UtLen = length(Ut);
   ApColor = L{StackCt,3};
   Xh=length(unique(ApColor));
   StackCt=StackCt-1;
     
   if length(Vt)+UtLen > CBCsize && Xh+UtLen > CBCsize
      ColorP=colorVertices(G(Vt,Vt));    % ColorP is current node's coloring; 
                                        % will be parent of nodes to be created next
      if UtLen > CBCsize
         CBCsize=UtLen;
         CBC = Ut;
      end
        
      if length(unique(ColorP)) + UtLen > CBCsize
         for i = Vt
            I1 = find(G(i,Vt)==0); %Indices of Vt that are neightbors of i
            N  = Vt(I1);           %Neighbors of i
            I2 = find(N<i);        %Indices of those neighbors with subscripts less than i
            Vtmp= N(I2);
            if ~isempty(Vtmp)
               S=sum(G(Vtmp,Vtmp));
               F=find(S==0);
               if ~isempty(F)
                  Utmp= [Ut i Vtmp(F)];
                  Vtmp=setdiff(Vtmp,Vtmp(F));
               else
                  Utmp=[Ut i];
               end
            else
               Utmp= [Ut i];
            end
               Xhtmp = ColorP(I1(I2));
               mXhtmp=max(Xhtmp);
            if isempty(mXhtmp)
               mXhtmp=0;
            end
            if mXhtmp + length(Utmp) > CBCsize && length(Vtmp) + length(Utmp) > CBCsize
               StackCt=StackCt+1;
               L{StackCt,1} = Vtmp; %#ok<AGROW>
               L{StackCt,2} = Utmp; %#ok<AGROW>
               L{StackCt,3} = Xhtmp;%#ok<AGROW>
            end 
         end 
      end 
   end 
end
  CBC=R(CBC);