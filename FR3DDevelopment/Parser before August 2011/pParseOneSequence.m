% pParseOneSequence(Node,Sequence) parses the single Sequence using the
% model in Node

% This program should only be called from pParseSequences, as that program
% initializes many required fields.

function [Sequence] = pParseOneSequence(Node,Sequence)

N = length(Node);         % number of nodes
x = Sequence.X;           % sequence to parse
L = length(x);            % length of sequence to parse

P = Sequence.P;

S(double('A')) = 1;
S(double('C')) = 2;
S(double('G')) = 3;
S(double('U')) = 4;
S(double('.')) = 5;

Code = S(x);                              % sequence of codes, not letters

t = cputime;

% ----------------------------------------- Omit certain nodes if desired
% Certain nodes can have length restrictions.  For instance, a
% hairpin node can only generate subsequences of length <= 5.  Nodes next to
% a hairpin can thus only generate subsequences of length <= 10, say.  There
% is no need to calculate the probability with which such a node generates
% much longer subsequences.


if isfield(Node(1),'minlength'),         % use length restrictions if desired
    picks=zeros(N,L);
    for n=1:N,
        if ~isempty(Node(n).minlength),
          Compute(n,:) = ((1:L) >= Node(n).minlength) .* ((1:L) <= Node(n).maxlength);
        else
          Compute(n,:) = ones(1,L);
        end
    end
else
    Compute = ones(N,L);
end
for p=1:L,
    NodesToCompute{p} = find(Compute(:,p));
end

% ----------------------------------------- note letter pair made by x(i),x(j)

for p=1:L,                                % length of subsequence to analyze
 for i=1:min(L,L-p+1),                    % starting position
  j=i+p-1;                                % ending position
  
  pc = Code(i) + (Code(j)-1)*4;           % paircode for [x(i) x(j)]
  if (Code(i) == 5) | (Code(j) == 5),     % for imaginary hairpins
    pc = 18;
  end
  for n=1:N,
      switch Node(n).type,
        case {'Initial','Hairpin','Junction','JunctionCluster'}
          P(n,1).transcode(i,j) = 1;     % must be active
        case 'Basepair'
          P(n,1).transcode(i,j) = 17;    % deleted
          P(n,2).transcode(i,j) = pc;    % active, with this paircode
        case 'Alternative'
          P(n,1).transcode(i,j) = 1;
          P(n,1).Alt(i,j)       = 1;      
        case 'Cluster'
          P(n,1).transcode(i,j) = 1;    % deleted
          P(n,2).transcode(i,j) = 2;    % active
      end
  end
 end                                     % loop over i
end                                      % loop over p

% --------------------------------------- inside/CYK algorithm

for p=1:L,                                % length of subsequence to analyze
 for i=1:min(L,L-p+1),                    % starting position
  j=i+p-1;                                % ending position
  if Sequence.Computed(i,j) == 0,         % if subseq has not been seen before
    pc = Code(i) + (Code(j)-1)*4;         % paircode for [x(i) x(j)]
    if (Code(i) == 5) | (Code(j) == 5),   % for imaginary hairpins
      pc = 18;
    end
    for m = length(NodesToCompute{p}):-1:1, % some nodes may be omitted  
      n = NodesToCompute{p}(m);             % current node to parse

      ok=1;

%[n Sequence.Constraint(n).iMin Sequence.Constraint(n).iMax Sequence.Constraint(n).jMin Sequence.Constraint(n).jMax]


      if (Sequence.Constraint(n).iMin <= i) && (i <= Sequence.Constraint(n).iMax) && (Sequence.Constraint(n).jMin <= j) && (j <= Sequence.Constraint(n).jMax),

      switch Node(n).type,
        % =================================================================
        case 'Initial'


            %          P = pParseInitial(Node,P,x,n,1,i,j);
            %function [P] = zParseInitial(node,P,x,n,s,i,j)
             s = 1;

             K = [[0 0]; [1 0]; [0 1]; [2 0]; [1 1]; [0 2]; [3 0]; [2 1]; [1 2]; ...
                  [0 3]; [3 1]; [2 2]; [1 3]; [3 2]; [2 3]; [3 3]]; % numbers of insertions

             mp  = -Inf;                            % maximum log prob for this n,s pair
             nextnum = Node(n).nextnode;

             if strcmp(Node(nextnum).type,'Alternative')
                 nextnum=Node(nextnum).nextnode(P(nextnum,1).Alt(i,j));
             end
             
             next = Node(nextnum);                  % next node in the tree
             maxmat = [];                           % matrix to store maximum values
                           % maxmat(h,c) is max prob using K(h,:) insertions on left
                           % and right and using child state c

             h = 1;                                % insertion pair we are processing
             while (h <= length(K(:,1))) & (sum(K(h,:)) <= j-i-1),
               a = i + K(h,1);                       % left side of child sequence
               b = j - K(h,2);                       % right side of child sequence
                   for c = 1:next.numstates,         % need to take care of alt node
                     tc = P(nextnum,c).transcode(a,b);   % code for transition prob
                     maxmat(h,c) = P(nextnum,c).mp(a,b) + next.lPIns(tc);
                   end
               h = h + 1;
             end

             h = h - 1;

             if h > 0,
               maxmat = maxmat + ((Node(n).loglip(s,1+K(1:h,1))+Node(n).logrip(s,1+K(1:h,2)))') * ones(1,next.numstates);
               [mv,ml] = max(maxmat);         % maximum over insertion combinations

               [mp,mc] = max(mv);  % add log probabilities for child states

               P(n,s).sub(i,j) = i + K(ml(mc),1);  % start of subsequence for next node
               P(n,s).sub(j,i) = j - K(ml(mc),2);  % end of subsequence for next node

               P(n,s).next(i,j) = mc;           % child state (of next node)
             end

            P(n,s).mp(i,j) = mp;                         % store the maximum prob

        case 'Hairpin'

            % ========================== Parse hairpin ============================

            %           P = pParseHairpin(Node,P,x,n,1,i,j);
            %function [P] = zParseHairpin(node,P,x,n,s,i,j)

            % In the future, the probabilities should depend on x(i:j)!
            % They should also depend on the hairpin type
            s = 1;

            if (j-i <= 4),
              switch Node(n).subtype
               case 'GNRA'
                 if (j-i == 3),
                   if ((x(j) == 'A') | ((x(i) ~= 'G') & (x(j) == 'C'))) & ...
                     (x(i+2) == 'A' | x(i+2) == 'G'),
                       P(n,s).mp(i,j) = log(10/904);
                   else 
                       P(n,s).mp(i,j) = log(1/904);
                   end
                 else
                     P(n,s).mp(i,j)   = log(1/(4+16+64+256));
                 end
               case 'NR'
                 if (j-i == 1),
                   if ((x(j) == 'A') | (x(j) == 'G')),
                     P(n,s).mp(i,j) = log(0.25*0.45);
                   else
                     P(n,s).mp(i,j) = log(0.25*0.05);
                   end 
                 end
               case '....'
                 P(n,s).mp(i,j) = -1000;
                 if (j-i == 3),
                   if (x(i) == '.') & (x(i+1)=='.') & (x(i+2)=='.') & (x(j)=='.'),
                     P(n,s).mp(i,j) = 1;
                   end
                 end
               case 'XX'
                 if (j-i == 1),
                   P(n,s).mp(i,j) = 0;
                 elseif (j-i < 4),
                   P(n,s).mp(i,j) = -3;
                 else
                   P(n,s).mp(i,j) = -10;
                 end
               case 'XXX'
                 if (j-i == 2),
                   P(n,s).mp(i,j) = 0;
                 elseif (j-i < 4),
                   P(n,s).mp(i,j) = -3;
                 else
                   P(n,s).mp(i,j) = -10;
                 end
               otherwise
                 P(n,s).mp(i,j)   = log(1/(4+16+64+256));
               end
            end

            P(n,s).sub(i,j)  = i;                     % start of subseq for next node
            P(n,s).sub(j,i)  = j;                     % end of subsequence
            P(n,s).next(i,j) = 0;                     % child state
            
        case 'Alternative'
            
            % ========================== Parse Alternative ============================
            nextnode = Node(n).nextnode;
            numAlt=length(nextnode);
            AltToUse=1;
            mp=-inf;
            state=1;

            for c=1:numAlt
                switch Node(nextnode(c)).type,
                   case {'Basepair'}
                       ss=2;
                   otherwise
                       ss=1;
                end
                for s=1:ss
                    m=P(nextnode(c),s).mp(i,j);
                    if m>mp
                        mp=m;
                        AltToUse=c;
                        state=s;
                    end
                end
            end
            P(n,1).mp(i,j)   = mp;               % store the maximum probability
            P(n,1).next(i,j) = state;
            P(n,1).sub(i,j)  = i;                  % start of subsequence for next node
            P(n,1).sub(j,i)  = j;                  % end of subsequence for next node           
            P(n,1).Alt(i,j)  = AltToUse;               % which Alternative to use
%[i j state AltToUse]


        case 'Basepair'
          if pc < 18 & i <= j-2,

            % ====================== Parse basepair ===================================

            nextnum = Node(n).nextnode;            % number of next node

            if strcmp(Node(nextnum).type,'Alternative')
                nextnum=Node(nextnum).nextnode(P(nextnum,1).Alt(i,j));
            end

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

            % ============================= End parse basepair

           end
        case 'Junction'

            % ============================= Parse junction
            % ----------------------------------- compute for node(n).state(s)
            % x(i,j) is divided as [i i+1 ... a-1 a a+1 ... j]
            %                       loop 1        loop 2

            nextnode = Node(n).nextnode;             % next node in the tree
          % need to add code here to compinsate for the alternative node
            for a=(i+1):j,                           % location of split 
              maxprob(a-i) = P(nextnode(1),1).mp(i,a-1) + P(nextnode(2),1).mp(a,j);
            end

            % also include the possibilities that each loop is deleted entirely
            % in the future

            [mv,ma] = max(maxprob);                  % most likely split

            P(n,s).sub(i,j)  = ma + i;               % location of the split
            P(n,s).next(i,j) = 0;                    % child state; not used
            P(n,s).mp(i,j)   = mv;                   % store the maximum prob
%[ma+i mv]
%pause
            % ----------------------------------- end computing for node(n).state(s)
        case 'JunctionCluster'

% need to add code here to compensate for the alternative node
            nextnum = Node(n).nextnode;            % number of next node
            next1   = Node(nextnum(1));               % next node in the tree
            next2   = Node(nextnum(2));               % next node in the tree
            maxmat  = [];                           % to store max values for insertions
            mp = -Inf;
            NL = length(Node(n).Left(1,:));       % number of bases on the left
            NR = length(Node(n).Right(1,:));
            NM = length(Node(n).Middle(1,:));
            aa = i;                             % in case the motif doesn't fit
            bb = j;
            cs = 1;
            lli = 1;
            rri = 1;
            cc  = i+1;
            mmi = 1;  
%             nml1=Node(n).minl(1);
%             nml2=Node(n).minl(2);
            mm  = i+1;
            for li = 1:length(Node(n).Left(:,1)),        % left insertion possibilities
              a = i + Node(n).Left(li,NL);               % first element for child
              for ri = 1:length(Node(n).Right(:,1)),     % right insertion possibilities
               b = j - Node(n).Right(ri,1);              % last element for child
               for m=(i+1):(j-1),                        % location of split between branches

                if (Sequence.Constraint(n).kMin <= m) && (m <= Sequence.Constraint(n).kMax),
                 for mi = 1:length(Node(n).Middle(:,1)),
                   c = m + Node(n).Middle(mi,NM);
                    if a <= m-1 & c<=b,
                    lc = Code(i - 1 + Node(n).Left(li,:));   % codes of bases on the left
                    mc = Code(m - 1 + Node(n).Middle(mi,:));  % codes of bases on the right
                    rc = Code(j + 1 - Node(n).Right(ri,:));  % codes of bases on the right
                    co = [lc mc rc];                            % codes of all bases to consider
                    if max(co) < 5,                          % letters include . for "hairpin"
                      S = Node(n).LLIP(li) + Node(n).LMIP(mi) + Node(n).LRIP(ri); % total score of interactions
                      for k = 1:length(Node(n).IBases(:,1)),   % loop through interactions
                        S = S+Node(n).Score(co(Node(n).IBases(k,1)),co(Node(n).IBases(k,2)),k);
                      end % end for k
                      maxprob= P(nextnum(1),1).mp(a,m-1) + P(nextnum(2),1).mp(c,b);
                      S = S + maxprob;
                      if S > mp,
                        mp  = S;
                        aa  = a;
                        bb  = b;
                        cc  = c;
                        lli = li;
                        rri = ri;
                        mmi = mi; 
                        mm  = m-1;
                      end %end if S > mp,
                    end % end if max(co) < 5, 
                   end % end if a <= b,
                 end % end for mi
                end

               end %for ri
              end % end for li
            end % end for m

            P(n,1).next(i,j) = 0;               % child state of next node
            P(n,1).next(j,i) = 0;               % child state of next node
            P(n,1).sub(i,j)  = aa;               % start of subsequence for next node
            P(n,1).sub(j,i)  = bb;               % end of subsequence for next node
            P(n,1).mp(i,j)   = mp;               % store the maximum probability
            P(n,1).mp(j,i)   = lli;          
            P(n,1).rmi(i,j)  = rri;
            P(n,1).rmi(j,i)  = mmi;
            P(n,1).sub2(i,j) = cc;
            P(n,1).sub2(j,i) = mm;

            % ------------------------------------ 
            
        case 'Cluster'

            % ============================= Parse cluster

            % zParseCluster(Node,P,Code,n,i,j) computes the probability that node n and its
            % children generate the subsequence x(i:j).  Information from lower nodes,
            % stored in the variable P, is required.

            nextnum = Node(n).nextnode;            % number of next node

            if strcmp(Node(nextnum).type,'Alternative')
                nextnum=Node(nextnum).nextnode(P(nextnum,1).Alt(i,j));
            end

            next    = Node(nextnum);               % next node in the tree

            maxmat = [];                           % to store max values for insertions

            % ------------------------------------ Compute probabilities if we delete

            for c = 1:next.numstates,              % loop through states of next node
              tc = P(nextnum,c).transcode(i,j);    % code for transition prob
              maxmat(1,c) = P(nextnum,c).mp(i,j) + next.lP(1,tc);
            end

            [mp,cs] = max(maxmat);                 % find the max prob child state

            P(n,1).mp(i,j)   = mp;                 % store the maximum probability
            P(n,1).next(i,j) = cs;                 % child state of next node
            P(n,1).sub(i,j)  = i;                  % start of subsequence for next node
            P(n,1).sub(j,i)  = j;                  % end of subsequence for next node

            % ------------------------------------ 

            mp = -Inf;
            maxmat = [];

            NL = length(Node(n).Left(1,:));              % number of bases on the left
            NR = length(Node(n).Right(1,:));

            aa = i;                                      % in case the motif doesn't fit
            bb = j;
            cs = 1;
            lli = 1;
            rri = 1;

            for li = 1:length(Node(n).Left(:,1)),        % left insertion possibilities
              a = i + Node(n).Left(li,NL);               % first element for child
              for ri = 1:length(Node(n).Right(:,1)),     % right insertion possibilities
               b = j - Node(n).Right(ri,1);              % last element for child
               if a <= b,                                % room for child too
                lc = Code(i - 1 + Node(n).Left(li,:));   % codes of bases on the left
                rc = Code(j + 1 - Node(n).Right(ri,:));  % codes of bases on the right
                co = [lc rc];                            % codes of all bases to consider
                if max(co) < 5,                          % letters include . for "hairpin"
                  S = Node(n).LLIP(li) + Node(n).LRIP(ri); % total score of interactions
                  for k = 1:length(Node(n).IBases(:,1)),   % loop through interactions
%[n li ri S k co(Node(n).IBases(k,1)) co(Node(n).IBases(k,2))]
                    S = S+Node(n).Score(co(Node(n).IBases(k,1)),co(Node(n).IBases(k,2)),k);
                  end

                  for c = 1:next.numstates,             % loop through states of next node
                    tc = P(nextnum,c).transcode(a,b);   % code for transition prob
                    maxmat(1,c) = P(nextnum,c).mp(a,b) + next.lPIns(tc);
                                                        % prob with this specific insertion
                                                        % no interaction between this motif
                                                        % and the next node
                  end

                  [y,mc] = max(maxmat);
                  S = S + y;

                  if S > mp,
                    mp  = S;
                    cs  = mc;
                    aa  = a;
                    bb  = b;
                    lli = li;
                    rri = ri;
                  end 
                else
                  S = [];
                end
               end
              end
            end

            P(n,2).mp(i,j)   = mp;               % store the maximum probability
            P(n,2).next(i,j) = cs;               % child state of next node
            P(n,2).sub(i,j)  = aa;               % start of subsequence for next node
            P(n,2).sub(j,i)  = bb;               % end of subsequence for next node

            P(n,2).mp(j,i)   = lli;
            P(n,2).next(j,i) = rri;              % store these somewhere; bad practice :(
            % ------------------------------------ 
      end
      end
    end
  end                                    % if not already computed
 end                                     % loop over i

% fprintf('Subsequences of length %3d took %5.2f minutes to analyze\n', p, (cputime-t)/60);
% drawnow
 t = cputime;

end                                      % loop over p

Sequence.P = P;

save('Parser/tempsave','Node','Sequence','L');
disp('------------------------')

Sequence.TraceInfo = pGetTraceInfo(Node,Sequence.P,L);

