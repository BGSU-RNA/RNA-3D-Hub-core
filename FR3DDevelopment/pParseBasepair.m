% zParseBasepair(node,P,Code,n,pc,i,j) computes the probability that node n 
% generates the subsequence x(i:j).  Information from lower nodes,
% stored in the variable P, is required.  
% The variable pc is the paircode for [x(i) x(j)].  It is passed in.
% Basepair nodes can be in two possible states:
% 1 - deleted - this node explains nothing, relying on the next node to
%               explain all of x(i:j)
% 2 - active  - this node generates x(i) and x(j) plus possibly some insertions

% The various combinations of the number of possible
% insertions on the right or left is the biggest user of time here.

% Later: Use the values in Code to determine the insertion probabilities!

function [P] = zParseBasepair(Node,P,Code,n,pc,i,j)

nextnum = Node(n).nextnode;            % number of next node
next    = Node(nextnum);               % next node in the tree

maxmat = [];                           % to store max values for insertions

% ------------------------------------ Compute probabilities if we delete

for c = 1:next.numstates,              % loop through states of next node
  tc = P(nextnum,c).transcode(i,j);    % code for transition prob
  maxmat(1,c) = P(nextnum,c).mp(i,j) + next.lP(17,tc);
end

[mp,cs] = max(maxmat);                 % find the max prob child state

P(n,1).mp(i,j)   = mp;                 % store the maximum probability
P(n,1).next(i,j) = cs;                 % child state of next node
P(n,1).sub(i,j)  = i;                  % start of subsequence for next node
P(n,1).sub(j,i)  = j;                  % end of subsequence for next node

% ------------------------------------- This state emits x(i) and x(j), no ins 

for c = 1:next.numstates,               % loop through states of next node
  tc = P(nextnum,c).transcode(i+1,j-1); % code for transition prob
  maxmat(1,c) = P(nextnum,c).mp(i+1,j-1) + next.lP(pc,tc);
end

%[c nextnum c i j pc tc]

[mp,cs] = max(maxmat);                  % find the max prob child state
mp = mp + Node(n).loglip(pc,1) + Node(n).logrip(pc,1); % add 0 insertion prob
aa = i + 1;                             % start of next node's sequence
bb = j - 1;                             % end of next node's sequence
 
% ------------------------------------- Analyze with insertions

K = [[1 0]; [0 1]; [2 0]; [1 1]; [0 2]; [3 0]; [2 1]; [1 2]; ...
     [0 3]; [3 1]; [2 2]; [1 3]; [3 2]; [2 3]; [3 3]]; 
                                        % numbers of insertions on L and R
                                        % sorted by total number (L+R)

% in the future: allow more insertions if the mean number of ins is high

h = 1;                                  % insertion pair we are processing
while (h <= length(K(:,1))) & (sum(K(h,:)) <= j-i-3),
  a = i + 1 + K(h,1);                   % left side of child sequence
  b = j - 1 - K(h,2);                   % right side of child sequence

  for c = 1:next.numstates,             % loop through states of next node
    tc = P(nextnum,c).transcode(a,b);   % code for transition prob
    maxmat(h,c) = P(nextnum,c).mp(a,b) + next.lPIns(tc);
                                        % prob with this specific insertion
  end

  h = h + 1;                            % go to next insertion combination
end

h = h - 1;                              % h is number of insertions analyzed

if (h > 0),                             % if any insertions fit in i:j
  maxmat = maxmat + ((Node(n).loglip(pc,1+K(1:h,1))+ ...
           Node(n).logrip(pc,1+K(1:h,2)))') * ones(1,next.numstates);
                                        % add probabilities of each number
                                        % of insertions on L and R

  [mv,ml] = max(maxmat);                % maximum over insertion combinations
                                        % gives a row of probabilities

  [nmp,mc] = max(mv);                   % maximum over child states

  if nmp > mp,                       % an insertion combination has higher prob
    mp = nmp;                           % new maximum probability
    cs = mc;                            % max prob child state
    aa = i + 1 + K(ml(mc),1);           % start of subsequence for next node
    bb = j - 1 - K(ml(mc),2);           % end of subsequence for next node
  end
end

P(n,2).mp(i,j)   = mp;               % store the maximum probability
P(n,2).next(i,j) = cs;               % child state of next node
P(n,2).sub(i,j)  = aa;               % start of subsequence for next node
P(n,2).sub(j,i)  = bb;               % end of subsequence for next node
