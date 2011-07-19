% pMakeConstraints(Node,Sequence) determines correspondences between the
% original 3D structure and the columns in the fasta file
% It also puts in substitution probabilities based on the observed interaction

function [Node,Sequence]=pMakeConstraints(Node,Sequence)

Minus=0;
N=length(Node);
Alt=0;
count=0;
Z=zeros(1,N);
ST=cputime;
NumAlt=1;
T=0;

for i=N:-1:1
  switch Node(i).type,
    case 'Initial'
%     Node(i).minlength=Node(i+1).minlength;
%     Node(i).maxlength=Node(i+1).maxlength+ceil(max(Node(i).lpar))+ceil(max(Node(i).rpar));
      L=Node(i).lpar(1,1);
      R=Node(i).rpar(1,1);
      X=0:10;
      Node(i).Z=sum(L.^X*exp(-L)./factorial(X)) * sum(R.^X*exp(-R)./factorial(X));  

    case 'Hairpin'

%     Node(i).minlength=1;
%     Node(i).maxlength=6;
      Z(i)=1;
    case 'Basepair'
      n1=pLtoN(Node(i).LeftLetter);
      n2=pLtoN(Node(i).RightLetter);
      Score        = pIsoScore(Node(i).Edge,n1,n2,Node(i).Delete);
      Node(i).P    = ones(17,1)*Score;         % state to state transitions
      Node(i).PIns = Score;                    % when no previous state           
  
      Minus=Minus+10*Node(i).P(17);
%     Node(i).minlength=Node(i+1).minlength + 2 - floor(Minus);
%     Node(i).maxlength=Node(i+1).maxlength + 1 + ceil(max(Node(i).lpar)) + ceil(max(Node(i).rpar));
      Minus=Minus-floor(Minus);
          
      L=Node(i).lpar(1,1);
      R=Node(i).rpar(1,1);
      X=0:10;     
      Node(i).Z=sum(L.^X*exp(-L)./factorial(X)) * sum(R.^X*exp(-R)./factorial(X));
          
    case 'Junction'
      [a,b]=Node(i).nextnode; 
%     Node(i).minlength=Node(a).minlength+Node(b).minlength;
%     Node(i).maxlength=Node(a).maxlength+Node(b).maxlength;
      Z(i)=1;



    case 'JunctionMotif'  % i will need to modify this
%     tic
      A=max(max(Node(i).IBases));
      V=ones(1,A);
      V(1)=0;
      for v=1:length(Node(i).Edge)
         n1=pLtoN(Node(i).Letters(Node(i).IBases(v,1)));
         n2=pLtoN(Node(i).Letters(Node(i).IBases(v,2)));
         Node(i).Score(:,:,v) = pIsoScore(Node(i).Edge(v),n1,n2);
      end

      IBases=Node(i).IBases;
      M=length(Node(i).Edge);
      Score=Node(i).Score;
      mm=4^A;
      Z(i)=0;
      for a=1:mm
              %----pCount----
              t=1;
                ok=0;
                while 1-ok
                    V(t)=V(t)+1;
                    ok=1;
                    if V(t)>4
                        V(t)=1;
                        t=t+1;
                        ok=0;
                    end
                end
%               V
              %----pCount----
              W=1;
              for v=1:M
                 W=W*exp(Score(V(IBases(v,1)),V(IBases(v,2)),v));
              end
              Z(i)=Z(i)+W;
          end
          Node(i).Score=Node(i).Score-log(Z(i))/M;
%           V
%           -log(Z(i))
%           toc
          a=Node(i).nextnode(1);  
          b=Node(i).nextnode(2);
%           Node(i).minlength=Node(a).minlength+Node(b).minlength;
%           Node(i).minl(1)=Node(b).minlength;
%           Node(i).minl(2)=Node(b).minlength;
%           Node(i).maxlength=Node(a).maxlength+Node(b).maxlength;        



    case 'Motif'
%           Node(i).minlength=Node(i+1).minlength+1;
%           Node(i).maxlength=Node(i+1).maxlength+1+max(max(Node(i).Left))+max(max(Node(i).Right));  
          A=max(max(Node(i).IBases));
          B=length(Node(i).IBases(:,1));
          IBases=Node(i).IBases;
          Score=Node(i).Score;
          V=ones(1,A);
          V(1)=0;
          mm=4^A;
          Z(i)=0;
          for a=1:mm
              %----pCount----
              t=1;
                ok=0;
                while 1-ok
                    V(t)=V(t)+1;
                    ok=1;
                    if V(t)>4
                        V(t)=1;
                        t=t+1;
                        ok=0;
                    end
                end
%               V=pCount(V); 
              %----pCount----
              W=1;
              for b=1:B
                  W=W*exp(Score(V(IBases(b,1)),V(IBases(b,2)),b));
              end 
              Z(i)=Z(i)+W;
          end
          Node(i).Score=Node(i).Score-log(Z(i))/B;




    case 'Alternative'
            % this is not quite right..fix it later
%           Node(i).minlength=min(cat(1,Node(Node(i).nextnode).minlength));
%           Node(i).maxlength=max(cat(1,Node(Node(i).nextnode).maxlength))+1; 
%           T=T+1;
%           NumAlt(T)=length(Node(i).nextnode);
          Z(i)=1;
  end


      if isfield(Node,'Alt')
          if Node(i).Alt > 0;
              count=count+1;
              Node(i).Alt;
            if Node(i).Alt==1
    %             [Alt i]
                for n=i:(m)
%                    Node(n).minlength=Node(n).minlength+Node(Alt+1).minlength-Node(m+1).minlength;
%                    Node(n).maxlength=Node(n).maxlength+Node(Alt+1).maxlength-Node(m+1).maxlength;
                end
            end
            if count==1
                Alt=i;
            elseif count==2
               m=i; 
            end
          end
      end

      % remove this later--------------------------
%                 Node(i).minlength=0;
%           Node(i).maxlength=122;
      %----------------------------------------------
      
      
end

% cputime-ST
% Z
% pause

% -log(sum(Z))
% loc=find(cat(1,Node.minlength)<1);
% for i=1:length(loc)
%     Node(loc(i)).minlength=1;
% end

if isfield(Node(1),'Bl')
    Sequence=SeqConstraints(Node,Sequence);
end


%--------------------------------------------------------------------
function Sequence=SeqConstraints(Node,Sequence)

M=length(Sequence);
N=length(Node);

F=zeros(M,N,3);
for k=1:M;
    
        L=length(Sequence(k).X);

        for n=1:N
%             [k n]
           if ~isempty(Node(n).Bl)
              q=0;
              w=length(Node(n).Bl);
              count=1;
              while q==0 & count<=w
                  [aa,f]=min(abs(Sequence(k).FastaIndex-Sequence(1).FastaIndex(Node(n).Bl(count))));

%                   f=find(Sequence(k).FastaIndex==Sequence(1).FastaIndex(Node(n).Bl(count)));
                  if length(f)>0
                      q=1;
                      Sequence(k).Constraint(n).Bl=Node(n).Bl-min(Node(n).Bl)+f;
                  end
                  count=count+1;
              end
              if q==0
                  if n==1
                      f=1;
                      Sequence(k).Constraint(n).Bl=Node(n).Bl-min(Node(n).Bl)+f;
                  else
                      f=F(k,n-1,1);
                      Sequence(k).Constraint(n).Bl=Node(n).Bl-min(Node(n).Bl)+f;
                  end

              end
              F(k,n,1)=f;    
           end


    %-----------------------
           if ~isempty(Node(n).Bm)
              q=0;
              w=length(Node(n).Bm);
              count=1;
              while q==0 & count<=w
                  [aa,f]=min(abs(Sequence(k).FastaIndex-Sequence(1).FastaIndex(Node(n).Bm(count))));
%                   f=find(Sequence(k).FastaIndex==Sequence(1).FastaIndex(Node(n).Bm(count)));
                  if length(f)>0
                      q=1;
                      Sequence(k).Constraint(n).Bm=Node(n).Bm-min(Node(n).Bm)+f;
                  end
                  count=count+1;
              end
              if q==0
                  f=F(k-1,n,2);
                  Sequence(k).Constraint(n).Bm=Node(n).Bm-min(Node(n).Bm)+f;
              end
              F(k,n,2)=f;   
           end

     %-------------------------
            if ~isempty(Node(n).Br)
              q=0;
              w=length(Node(n).Br);
              count=1;
              while q==0 & count<=w
                  [aa,f]=min(abs(Sequence(k).FastaIndex-Sequence(1).FastaIndex(Node(n).Br(count))));
%                   f=find(Sequence(k).FastaIndex==Sequence(1).FastaIndex(N
%                   ode(n).Br(count)));
                  if length(f)>0
                      q=1;
                      Sequence(k).Constraint(n).Br=Node(n).Br-min(Node(n).Br)+f;
                  end
                  count=count+1;
              end
              if q==0
                  if n==1
                      f=L;
                      Sequence(k).Constraint(n).Br=Node(n).Br-min(Node(n).Br)+f;
                  else
                      f=F(k,n-1,3);
                      Sequence(k).Constraint(n).Br=Node(n).Br-min(Node(n).Br)+f;
                  end
              end
              F(k,n,3)=f;    
           end
        end
 
        
            
    Sequence(k).Constraint=pConstructBases(Node,Sequence(k).Constraint,L);

    picks=zeros(N,L);
    for n=1:N,
        if ~isempty(Sequence(k).Constraint(n).Bases),
          Bases=Sequence(k).Constraint(n).Bases(find(Sequence(k).Constraint(n).Bases<=L));
          picks(n,:)=zeros(1,L);
          picks(n,Bases)=ones(1,length(Bases));
        else
          picks(n,:)=ones(1,L);
        end
    end
    Sequence(k).ToCompute=sparse(picks);
end
    
%--------------------------------------------------------------------
% function V=pCount(V)
% 
% t=1;
% ok=0;
% while 1-ok
%     V(t)=V(t)+1;
%     ok=1;
%     if V(t)>4
%         V(t)=1;
%         t=t+1;
%         ok=0;
%     end
% end
