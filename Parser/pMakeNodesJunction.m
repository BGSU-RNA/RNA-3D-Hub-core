% pMakeNodesJunction identifies one or more junctions and make nodes for them

      C1 = full(sum(sum(H(r:r+jcdepth,t:t+jcdepth))));   % junction cluster 1-3
      C2 = full(sum(sum(H(r:r+jcdepth:u-jcdepth:u))));   % junction cluster 1-4
      C3 = full(sum(sum(H(s-jcdepth:s,t:t+jcdepth))));   % junction cluster 2-3
      C4 = full(sum(sum(H(s-jcdepth:s,u-jcdepth:u))));   % junction cluster 2-4

if C1+C2+C3+C4 > 0,
% [C1 C2 C3 C4]
  disp('Junction includes some nested pairs that will be removed.');
end



      C1 = sum(sum(G(r:r+jcdepth,t:t+jcdepth)~=0));   % junction cluster 1-3
      C2 = sum(sum(G(r:r+jcdepth:u-jcdepth:u)~=0));   % junction cluster 1-4
      C3 = sum(sum(G(s-jcdepth:s,t:t+jcdepth)~=0));   % junction cluster 2-3
      C4 = sum(sum(G(s-jcdepth:s,u-jcdepth:u)~=0));   % junction cluster 2-4

      [i,j,k] = find(G(r:r+jcdepth,t:t+jcdepth));   % junction cluster 1-3
      kk = k;
      [i,j,k] = find(G(r:r+jcdepth:u-jcdepth:u));   % junction cluster 1-4
      kk = [kk; k];
      [i,j,k] = find(G(s-jcdepth:s,t:t+jcdepth));   % junction cluster 2-3
      kk = [kk; k];
      [i,j,k] = find(G(s-jcdepth:s,u-jcdepth:u));   % junction cluster 2-4
      kk = [kk; k];

Node(1).JunctionDeletion = [Node(1).JunctionDeletion; kk];

if Verbose > 1,
  full([C1 C2 C3 C4])
end

      % ------------------ Remove junction clusters --------------------
      % Note:  this is only happening between the two loops identified so
      % far, but there may be more loops, and between them these interactions
      % are not being removed!  Maybe this is not actually a problem.
      % The little study above tells me that in E. coli 16S, no nested pairs
      % were removed by this.

      if C1 > 0,
        G(r:r+jcdepth,t:t+jcdepth) = 0*G(r:r+jcdepth,t:t+jcdepth);
      end
      if C2 > 0,
        G(r:r+jcdepth:u-jcdepth:u) = 0*G(r:r+jcdepth:u-jcdepth:u);
      end
      if C3 > 0,
        G(s-jcdepth:s,t:t+jcdepth) = 0*G(s-jcdepth:s,t:t+jcdepth);
      end
      if C4 > 0,
        G(s-jcdepth:s,u-jcdepth:u) = 0*G(s-jcdepth:s,u-jcdepth:u);
      end

      C1 = full(sum(sum(G(r:r+jcdepth,t:t+jcdepth))));   % junction cluster 1-3
      C2 = full(sum(sum(G(r:r+jcdepth:u-jcdepth:u))));   % junction cluster 1-4
      C3 = full(sum(sum(G(s-jcdepth:s,t:t+jcdepth))));   % junction cluster 2-3
      C4 = full(sum(sum(G(s-jcdepth:s,u-jcdepth:u))));   % junction cluster 2-4

      if C1+C2+C3+C4 == 0,                  % no interaction across junction

        junc = [];                          % indices where loops start & end

        while (sum(sum(H((r+1):(s-1),(r+1):(s-1)) > 0)) > 0) && ...
              (sum(sum(H((s+1):(u),(s+1):(u)) > 0)) > 0),  % still two loops

          % probe for start of next block of nested pairs, 
          % s+1 to t-1 is junction strand
          while sum(H(t,(t+1):B) > 0) == 0 && t < u,
                                         % if t does not make a nested pair,
            t = t + 1;
          end

          % decide which loop gets the strand between the loops
          % these are the indices that it could interact with in a cluster

          rrr = unique([rr:min(r+cdepth,s) max(rr,s-cdepth):s]);
          ttt = unique([t:min(t+cdepth,u) max(Interact{t}.Index(1)-cdepth,t):u]);

          % here the criterion is simple, just choose whichever has the
          % larger number of interactions.  But it would be better yet to
          % split the loop, as with 556:567 in 2avy

          if sum(sum(G(rrr,(s+1):(t-1))>0)) >= sum(sum(G((s+1):(t-1),ttt))),
            junc = [junc; [rr t-1]];           % store limits of this loop
          else
            junc = [junc; [rr s]];
          end

          r  = junc(end,2) + 1;           % NT after the first loop ends
          rr = r;                         % copy of that, for probing forward
          s  = Interact{t}.Index(1);      % the far side of the next loop
          t  = s + 1;                     % one beyond that
        end

        junc = [junc; [r u]];             % store limits of this last loop

        NL = length(junc(:,1));             % number of loops

        id = fix(10000*rand);

        if Verbose > 0,

%[File.NT(a).Number ' ' File.NT(b).Number ' ' File.NT(B).Number]

          fprintf('\nJunction with %d loops, call it J%d\n', NL,id);
          for ln = 1:NL,
            fprintf('Actual loop %d of junction J%d - Nucleotides %5s_%s to %5s_%s, length %3d\n',ln,id,File.NT(junc(ln,1)).Number,File.NT(junc(ln,1)).Chain,File.NT(junc(ln,2)).Number,File.NT(junc(ln,2)).Chain,junc(ln,2)+1-junc(ln,1));
          end
        end
  
        n = n+1;                              % move to next node
        Node(n).type       = 'Junction';      % junction with no cluster
        Node(n).LeftIndex  = a;
        Node(n).RightIndex = B;
        Node(n).NumLoops   = 2;
        Node(n).id         = ['J' num2str(id)];

if Verbose > 0,
  fprintf('%3d Junction\n', n);
end

        jn = n;                               % index of this node

        if NL == 2,                           % exactly two branches,

         for ln = 1:NL,
          if Verbose > 0,
            fprintf('\n');
            fprintf('Loop %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',ln,NL,id,File.NT(junc(ln,1)).Number,File.NT(junc(ln,2)).Number,junc(ln,2)+1-junc(ln,1));
          end
          nn = length(Node) + 1;
          Node(jn).nextnode(ln) =  length(Node)+1;
          Node = pMakeNodes(File,Param,junc(ln,1),junc(ln,2),Truncate,Data,Node,n);
          Node(nn).id         = ['J' num2str(id)];
          n = length(Node);
         end


        else                                    % more than two branches

          NN = ceil(NL/2);                      % # branches for 1st child

          if Verbose > 0,
            fprintf('\n');
            fprintf('Loop %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',1,2,id,File.NT(junc(1,1)).Number,File.NT(junc(NN,2)).Number,junc(NN,2)+1-junc(1,1));
          end

          nn = length(Node) + 1;
          Node(jn).nextnode(1) = length(Node) + 1;
          Node = pMakeNodes(File,Param,junc(1,1),junc(NN,2),Truncate,Data,Node,n);
          Node(nn).id         = ['J' num2str(id)];
          n = length(Node);

          if Verbose > 0,
            fprintf('\n');
            fprintf('Loop %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',2,2,id,File.NT(junc(NN+1,1)).Number,File.NT(junc(NL,2)).Number,junc(NL,2)+1-junc(NN+1,1));
          end

          nn = length(Node) + 1;
          Node(jn).nextnode(2) = length(Node) + 1;
          Node = pMakeNodes(File,Param,junc(NN+1,1),junc(NL,2),Truncate,Data,Node,n);
          Node(nn).id         = ['J' num2str(id)];

        end

        Node(n).Comment = ['// Junction node ' File.NT(Node(n).LeftIndex).Base File.NT(Node(n).LeftIndex).Number ' - ' File.NT(Node(n).RightIndex).Base File.NT(Node(n).RightIndex).Number ' ID ' Node(n).id];

      else                                      % two-loop junction cluster

        % find extent of interactions between loops

        t = b;

        rr = r + jcdepth;
        while (sum(G(rr,[t:t+jcdepth u-jcdepth:u])) == 0) && (rr > r),
          rr = rr - 1;
        end

        ss = s - jcdepth;
        while (sum(G(ss,[t:t+jcdepth u-jcdepth:u])) == 0) && (ss < s),
          ss = ss + 1;
        end

        tt = t + jcdepth;
        while (sum(G([r:r+jcdepth s-jcdepth:s],tt)) == 0) && (tt > t),
          tt = tt - 1;
        end

        uu = u - jcdepth;
        while (sum(G([r:r+jcdepth s-jcdepth:s],uu)) == 0) && (uu < u),
          uu = uu + 1;
        end

        % second, extent of additional interactions within loops

        rrr = r + jcdepth;
        while (sum(G(rrr,[r:rr ss:s])) == 0) && (rrr > rr),
          rrr = rrr - 1;
        end

        sss = s - jcdepth;
        while (sum(G(sss,[r:rr ss:s])) == 0) && (sss < ss),
          sss = sss + 1;
        end

        ttt = t + jcdepth;
        while (sum(G([t:tt uu:u],ttt)) == 0) && (ttt > tt),
          ttt = ttt - 1;
        end

        uuu = u - jcdepth;
        while (sum(G([t:tt uu:u],uuu)) == 0) && (uuu < uu),
          uuu = uuu + 1;
        end

        n = n + 1;
        Node(n).type        = 'JunctionCluster';  % 
        Node(n).LeftIndex   = [r:rrr];
        Node(n).MiddleIndex = [sss:ttt];
        Node(n).RightIndex  = [uuu:u];

%        Node(n).Left(1,:)   = union(zs,xs);
%        Node(n).Middle(1,:) = fliplr(union(zt,yt)); % correct?
%        Node(n).Right(1,:)  = fliplr(union(zt,yt)); % correct?

        Node(n).LIP = [1];
        Node(n).MIP = [1];
        Node(n).RIP = [1];

        % add additional insertion combinations and probabilities here!
        % add scores for the various basepairs here!

        if Verbose > 0,
          fprintf('%3d Junction Cluster %4s %4s %4s %4s %4s %4s\n', n, File.NT(r).Number, File.NT(rrr).Number, File.NT(sss).Number, File.NT(ttt).Number, File.NT(uuu).Number, File.NT(u).Number);
          fprintf('================================================================================================================\n');
        end

        r = rrr;
        s = sss;
        t = ttt;
        u = uuu;

        Node(n).nextnode(1) =  n+1;          % index of next node in tree
        Node = pMakeNodes(File,Param,r,s,Truncate,Data,Node,n);

        Node(n).nextnode(2)  = length(Node)+1;
        Node = pMakeNodes(File,Param,t,u,Truncate,Data,Node,length(Node));
      end                                  % junction cluster

      Node(n).P    = [0.05*ones(17,1) 0.95*ones(17,1)];
                                            % state to state transitions
      Node(n).PIns = [0.05 0.95];   % when no previous state

      EndLoop = 1;

