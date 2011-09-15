% pParseJunctionMotif(Node,P,Code,n,i,j) computes the probability that node n and its
% children generate the subsequence x(i:j).  Information from lower nodes,
% stored in the variable P, is required.
function [P] = pParseJunctionMotif(Node,P,Code,n,i,j)
nextnum = Node(n).nextnode;            % number of next node
next1    = Node(nextnum(1));               % next node in the tree
next2    = Node(nextnum(2));               % next node in the tree
maxmat = [];                           % to store max values for insertions
% % ------------------------------------ Compute probabilities if we delete
% 
% for c = 1:next.numstates,              % loop through states of next node
%   tc = P(nextnum,c).transcode(i,j);    % code for transition prob
%   maxmat(1,c) = P(nextnum,c).mp(i,j) + next.lP(1,tc);
% end
% 
% [mp,cs] = max(maxmat);                 % find the max prob child state
% 
% P(n,1).mp(i,j)   = mp;                 % store the maximum probability
% P(n,1).next(i,j) = cs;                 % child state of next node
% P(n,1).sub(i,j)  = i;                  % start of subsequence for next node
% P(n,1).sub(j,i)  = j;                  % end of subsequence for next node
% 
% % ------------------------------------ 
mp = -Inf;
NL = length(Node(n).Left(1,:));              % number of bases on the left
NR = length(Node(n).Right(1,:));
aa = i;                                      % in case the motif doesn't fit
bb = j;
cs = 1;
lli = 1;
rri = 1;
%-----------------------------------------------------
for m=(i+1):(j-1),
%-----------------------------------------------------
for li = 1:length(Node(n).Left(:,1)),        % left insertion possibilities
  a = i + Node(n).Left(li,NL);               % first element for child
  for ri = 1:length(Node(n).Right(:,1)),     % right insertion possibilities
   b = j - Node(n).Right(ri,1);              % last element for child
   for mi = 1:length(Node(n).Middle(:,1),
       c = m + Node(n).Middle(mi,1);
       if a <= c & c<=b,                                % room for child too
        lc = Code(i - 1 + Node(n).Left(li,:));   % codes of bases on the left
        mc = Code(m - 1 + Node(n).Middle(mi,:));  % codes of bases on the right
        rc = Code(j + 1 - Node(n).Right(ri,:));  % codes of bases on the right
        co = [lc mc rc];                            % codes of all bases to consider
        if max(co) < 5,                          % letters include . for "hairpin"
          S = Node(n).LLIP(li) + Node(n).LMIP(mi) + Node(n).LRIP(ri); % total score of interactions
          for k = 1:length(Node(n).IBases(:,1)),   % loop through interactions
            S = S+Node(n).Score(co(Node(n).IBases(k,1)),co(Node(n).IBases(k,2)),k);
          end % end for k
          for d = 1:next.numstates,             % loop through states of next node
            td1 = P(nextnum,d).transcode(a,c);   % code for transition prob
            maxmat1(d) = P(nextnum(1),d).mp(a,c) + next.lPIns(td1);
                                                % prob with this specific insertion
                                                % no interaction between this motif
                                                % and the next node
            td2 = P(nextnum,d).transcode(c,b);   % code for transition prob
            maxmat2(d) = P(nextnum(2),d).mp(c,b) + next.lPIns(td2);
                                                % prob with this specific insertion
                                                % no interaction between this motif
                                                % and the next node                                                
          end %end for d
          [y1,mc1] = max(maxmat1);
          [y2,mc2] = max(maxmat2);
          S = S + y1 + y2;
          if S > mp,
            mp  = S;
            cs1 = mc1;
            cs2 = mc2;
            aa  = a;
            bb  = b;
            cc  = c;
            lli = li;
            rri = ri;
            P(n,2).mp(i,j)   = mp;               % store the maximum probability
            % P(n,2).next(i,j) = cs;               % child state of next node
            P(n,2).next(i,j) = 0;               % child state of next node
            % P(n,2).sub(i,j)  = aa;               % start of subsequence for next node
            % P(n,2).sub(j,i)  = bb;               % end of subsequence for next node
            P(n,2).sub(i,j)  = m;               % start of subsequence for next node
%             P(n,2).sub(j,i)  = bb;               % end of subsequence for next node
%             P(n,2).mp(j,i)   = lli;
%             P(n,2).next(j,i) = rri;              % store these somewhere; bad practice :( 
          end %end if S > mp,
        end % end if max(co) < 5, 
       end % end if a <= b,
   end % end for c
  end % end for ri
end % end for li
end %for m
% ------------------------------------ 