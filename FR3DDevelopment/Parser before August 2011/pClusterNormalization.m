% pClusterNormalization calculates the normalization factor for cluster probabilities

function [N] = pClusterNormalization(N)

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

Node(i).Score = Node(i).Score-log(Z(i))/B;  
