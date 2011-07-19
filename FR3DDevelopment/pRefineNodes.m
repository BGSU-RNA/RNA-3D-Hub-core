% pRefineNodes(Node,Sequence) determines correspondences between the
% original 3D structure and the columns in the fasta file
% It also puts in substitution probabilities based on the observed interaction

function [Node,Sequence] = pRefineNodes(Node,Sequence)

N = length(Node);

Alt=0;
count=0;
Z=zeros(1,N);
ST=cputime;
NumAlt=1;
T=0;

for i = length(Node):-1:1,                            % loop through nodes

  switch Node(i).type,
    case 'Initial'

    case 'Hairpin'

    case 'Basepair'
          
    case 'Junction'

    case 'JunctionCluster'  % i will need to modify this

      A = max(max(Node(i).IBases));
      V = ones(1,A);
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

    case 'Cluster'
      A = max(max(Node(i).IBases));      % size of cluster
      B = length(Node(i).IBases(:,1));   % number of interactions
      IBases=Node(i).IBases;
      Score=Node(i).Score;
      V=ones(1,A);
      V(1)=0;
      mm=4^A;                            % mm = 4 to number of bases
      Z(i)=0;                            % normalization constant
      for a=1:mm                         % run through all possibilities
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
%       V=pCount(V); 

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

  end                                             % switch Node(i).type


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
      
      
end                                              % loop through nodes

% cputime-ST
% Z
% pause

% -log(sum(Z))
% loc=find(cat(1,Node.minlength)<1);
% for i=1:length(loc)
%     Node(loc(i)).minlength=1;
% end

