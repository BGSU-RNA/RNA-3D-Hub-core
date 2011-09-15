% pGetTraceInfo(Node,P,L) extracts trace info from P

function [t] = pGetTraceInfo(Node,P,L)

N = length(Node);
Active=0;
t(1).i     = 1;
t(1).j     = L;                              % length of sequence
t(1).state = 1;
t(1).active= 1;
for k = 2:length(Node),
    t(k).active= 0;
end
for n = 1:length(Node),
%   [n Active]
  nextnode = Node(n).nextnode;
  i = t(n).i;
  j = t(n).j;
  s = t(n).state;
 if  t(n).active==1
  switch Node(n).type,
    case {'Initial','Basepair'}

%[n s i j t(nextnode).active Active]
          
      t(n).a            = P(n,s).sub(i,j);
      t(n).b            = P(n,s).sub(j,i);
      t(n).mp           = P(n,s).mp(i,j);
      if t(nextnode).active==0 & n~=Active
          t(nextnode).i     = t(n).a;
          t(nextnode).j     = t(n).b;
          t(nextnode).state = P(n,s).next(i,j);
          t(nextnode).active= max(t(n).active,t(nextnode).active);
      end
    case {'Cluster'}
%   [n s i j]

      t(n).a            = P(n,s).sub(i,j);
      t(n).b            = P(n,s).sub(j,i);
      t(n).l            = P(n,s).mp(j,i);
      t(n).r            = P(n,s).next(j,i);
      t(n).mp           = P(n,s).mp(i,j);
      if (t(nextnode).active==0) && (n ~= Active),
          t(nextnode).i     = t(n).a;
          t(nextnode).j     = t(n).b;
          t(nextnode).state = P(n,s).next(i,j);
          t(nextnode).active= max(t(n).active,t(nextnode).active);
      end
    case 'Hairpin'
%         [n s i j]
      t(n).mp           = P(n,s).mp(i,j);
    case 'Alternative'
%          [n s i j]
%         mp=P(n,1).mp(i,j);
%         mp
%         [n s i j]
%         P(n,1).Alt(i,j)
%         [n s i j mp P(n,1).Alt(i,j)]
%         [P(n,1).sub(i,j) P(n,1).sub(j,i)]

        t(n).mp               = P(n,1).mp(i,j);
        t(n).AltToUse         = P(n,1).Alt(i,j);
        t(n).a                = P(n,1).sub(i,j);
        t(n).b                = P(n,1).sub(j,i);
        for k=1:length(Node(n).nextnode)
           t(Node(n).nextnode(k)).state = P(n,1).next(i,j);
           t(Node(n).nextnode(k)).i     = t(n).a;
           t(Node(n).nextnode(k)).j     = t(n).b;     
        end
        
        t(nextnode(P(n,1).Alt(i,j))).active= 1;
        if t(n).AltToUse==length(Node(n).nextnode)
           Active=1;
        else
            Active=nextnode(t(n).AltToUse);
            while Active+1==Node(Active).nextnode(1)
             Active=Active+1;
            end
            Active=Node(Active).nextnode(1);
        end
        Active=Active-1;
 
    case 'Junction' 
%     [n s i j]
     if t(nextnode(1)).active==0 & n~=Active
      t(nextnode(1)).i     = t(n).i;
      t(nextnode(1)).j     = uint8(double(P(n,s).sub(i,j))-1);
      t(nextnode(1)).state = 1;
      t(nextnode(2)).i     = P(n,s).sub(i,j);
      t(nextnode(2)).j     = t(n).j;
      t(nextnode(2)).state = 1;
      t(nextnode(1)).active= max(t(n).active,t(nextnode(1)).active);
      t(nextnode(2)).active= max(t(n).active,t(nextnode(2)).active);
     end
    case 'JunctionCluster' 
%     [n s i j]   
      t(n).a               = P(n,s).sub(i,j);
      t(n).b               = P(n,s).sub(j,i);
%       t(n).m               = uint8(double(P(n,s).sub2(j,i))-1);
      t(n).m               = P(n,s).sub2(j,i);
      t(n).c               = P(n,s).sub2(i,j);
      t(n).l               = P(n,s).mp(j,i);
      t(n).r               = P(n,s).rmi(i,j);
      t(n).cm              = P(n,s).rmi(j,i);
      t(n).mp              = P(n,s).mp(i,j);
      if t(nextnode(1)).active==0 & n~=Active
          t(nextnode(1)).i     = t(n).a;
          t(nextnode(1)).j     = t(n).m;
          t(nextnode(2)).i     = t(n).c;
          t(nextnode(2)).j     = t(n).b;
          t(nextnode(1)).state = 1;                
          t(nextnode(2)).state = 1;
          t(nextnode(1)).active= max(t(n).active,t(nextnode(1)).active);
          t(nextnode(2)).active= max(t(n).active,t(nextnode(2)).active);
      end
  end
 end
end