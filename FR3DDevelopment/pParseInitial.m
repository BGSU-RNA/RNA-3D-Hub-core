function [P] = zParseInitial(node,P,x,n,s,i,j)

 K = [[0 0]; [1 0]; [0 1]; [2 0]; [1 1]; [0 2]; [3 0]; [2 1]; [1 2]; ...
      [0 3]; [3 1]; [2 2]; [1 3]; [3 2]; [2 3]; [3 3]]; % numbers of insertions

 mp  = -Inf;                            % maximum log prob for this n,s pair

 nextnum = node(n).nextnode;
 next = node(nextnum);                  % next node in the tree

 maxmat = [];                           % matrix to store maximum values
               % maxmat(h,c) is max prob using K(h,:) insertions on left
               % and right and using child state c

 h = 1;                                % insertion pair we are processing
 while (h <= length(K(:,1))) & (sum(K(h,:)) <= j-i-1),
   a = i + K(h,1);                       % left side of child sequence
   b = j - K(h,2);                       % right side of child sequence

   for c = 1:next.numstates,
     tc = P(nextnum,c).transcode(a,b);   % code for transition prob
     maxmat(h,c) = P(nextnum,c).mp(a,b) + next.lPIns(tc);
   end
   
   h = h + 1;
 end

 h = h - 1;

 if h > 0,
   maxmat = maxmat + ((node(n).loglip(s,1+K(1:h,1))+node(n).logrip(s,1+K(1:h,2)))') * ones(1,next.numstates);
   [mv,ml] = max(maxmat);         % maximum over insertion combinations

   [mp,mc] = max(mv);  % add log probabilities for child states

   P(n,s).sub(i,j) = i + K(ml(mc),1);  % start of subsequence for next node
   P(n,s).sub(j,i) = j - K(ml(mc),2);  % end of subsequence for next node

   P(n,s).next(i,j) = mc;           % child state (of next node)
 end

P(n,s).mp(i,j) = mp;                         % store the maximum prob

