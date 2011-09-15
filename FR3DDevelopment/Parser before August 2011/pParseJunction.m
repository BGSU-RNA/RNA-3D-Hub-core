% zParseJunction parses a junction.  It assumes that the child nodes are
% Initial nodes, and so only have one state they can be in.

function [P] = zParseJunction(node,P,x,n,s,i,j)

% ----------------------------------- compute for node(n).state(s)
% x(i,j) is divided as [i i+1 ... a-1 a a+1 ... j]
%                       loop 1        loop 2

nextnode = node(n).nextnode;             % next node in the tree

for a=(i+1):j,                           % location of split 
  maxprob(a-i) = P(nextnode(1),1).mp(i,a-1) + P(nextnode(2),1).mp(a,j);
end

% also include the possibilities that each loop is deleted entirely
% in the future

[mv,ma] = max(maxprob);                  % most likely split

P(n,s).sub(i,j)  = ma + i;               % location of the split
P(n,s).next(i,j) = 0;                    % child state; not used
P(n,s).mp(i,j)   = mv;                   % store the maximum prob

% ----------------------------------- end computing for node(n).state(s)
