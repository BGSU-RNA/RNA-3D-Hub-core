% pMakeNodes(File,NTNumber) makes a secondary structure node model based on
% the Edge interaction matrix in File, starting at NTNumber and ending at
% LastNTNumber.  It assigns various nodes consistent with this secondary
% structure.  

function [Node] = pMakeNodes(File,NTNumber,LastNTNumber,Interact,Node,n)

verbose = 1;

cdepth  = 7;                      % how far to look ahead for a cluster
jcdepth = 4;                      % how far to look for a junction cluster

% if File is a text string (filename), load the file and display

if nargin<6
   n=0;
end

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char'),
  NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell'),
  Indices = zIndexLookup(File,NTNumber);
else
  Indices = NTNumber;
end

N = length(File.NT);                       % number of nucleotides
E = abs(fix(File.Edge));                   % don't distinguish subcategories
G = E .* (E < 15) .* (E ~= 0);             % consider basepairing only

if nargin == 2,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
end

if nargin < 4,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
  [y,L] = sort(E(a,k));                    % sort by edge interaction category
  Interact{a}.Categ = abs(File.Edge(a,k(L)));   % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

a  = Indices(1);                           % first index; current index
A  = a;                                    % previous interacting base on left
AA = a;                                    % previous cWW base on left

B  = LastNTNumber;                         % next base on right
BB = a;                                    % previous cWW base on right

if verbose > 0,
  fprintf('Loop %4s %4s\n', File.NT(a).Number, File.NT(B).Number);
end

% Initial node creation ------

n = n+1;                                   % move to next node
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).lpar      = 0.01;                  % left insertion parameter
Node(n).rpar      = 0.01;         % ARE THESE STILL USED?
Node(n).LeftIndex = a;
Node(n).RightIndex= B;
L = Node(n).lpar(1,1);
R = Node(n).rpar(1,1);
X = 0:10;
Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
          * sum(R.^X*exp(-R)./factorial(X));%normalization constant

                                           % probe for insertions
while sum(G(a,a+1:B)) == 0,                % no interaction w/in this loop
  Node(n).lpar = Node(n).lpar + 1;         % increase mean number of insertions
  a = a + 1;                               % move to next base
end

while sum(G(a:B-1,B)) == 0,                % no interaction w/in this loop
  Node(n).rpar = Node(n).rpar + 1;         % increase mean number of insertions
  B = B - 1;                               % move to next base
end


Node(n).leftLengthDist = subPoisson(Node(n).lpar);
Node(n).leftLetterDist = [0.25 0.25 0.25 0.25];
Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

if verbose > 0,
  fprintf('%3d Initial %s:%s and %s:%s\n', n, File.NT(Node(n).LeftIndex).Number, File.NT(a).Number, File.NT(Node(n).RightIndex).Number, File.NT(B).Number);
end

% -----------------------------

EndLoop = 0;                               % flag for the end of the loop

while (EndLoop == 0) & (a <= N),           % while not the end of the loop,

  b = Interact{a}.Index(1);              % index of what a interacts with

  if (a < b),                            % if b comes after a

    % ---------------------------------- Check for junction

    % check to see if a and b are now in a new loop with cWW's

    r = a;
    s = b;
    t = b+1;
    u = B;

%[File.NT(a).Number ' ' File.NT(b).Number ' ' File.NT(B).Number]
%[a b B]


    if (sum(sum(G(t:u,t:u) == 1))         > 0) && ...
       (sum(sum(G(a+1:b-1,a+1:b-1) == 1)) > 0),
                    % there are helices between a and b and between b+1 and B

      C1 = full(sum(sum(G(r:r+jcdepth,t:t+jcdepth))));   % junction cluster 1-3
      C2 = full(sum(sum(G(r:r+jcdepth:u-jcdepth:u))));   % junction cluster 1-4
      C3 = full(sum(sum(G(s-jcdepth:s,t:t+jcdepth))));   % junction cluster 2-3
      C4 = full(sum(sum(G(s-jcdepth:s,u-jcdepth:u))));   % junction cluster 2-4

%[C1 C2 C3 C4]

      if C1+C2+C3+C4 == 0,                  % no interaction across junction

        junc = [];                          % indices where loops start & end

        while (sum(sum(G((r+1):(s-1),(r+1):(s-1)) == 1)) > 0) && ...
              (sum(sum(G((s+1):(u),(s+1):(u)) == 1)) > 0),

          while (sum(G(s,(s+1):(u))) == 0) && (s < u),
            s = s + 1;
          end
          junc = [junc; [r s-1]];           % store limits of this loop

          r = s;                            % move to next loop
          s = Interact{r}.Index(1);        
        end

        junc = [junc; [r u]];             % store limits of this last loop

        NL = length(junc(:,1));             % number of loops

        if verbose > 0,

%[File.NT(a).Number ' ' File.NT(b).Number ' ' File.NT(B).Number]

          id = fix(10000*rand);
          fprintf('\nJunction with %d loops, call it J%d\n', NL,id);
          for ln = 1:NL,
            fprintf('Loop %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',ln,id,File.NT(junc(ln,1)).Number,File.NT(junc(ln,2)).Number,junc(ln,2)+1-junc(ln,1));
          end
        end
  
        n = n+1;                              % move to next node
        Node(n).type       = 'Junction';       % junction with no cluster
        Node(n).LeftIndex  = a;
        Node(n).RightIndex = B;
        Node(n).NumLoops   = NL;

        jn = n;                               % index of this node
%n
        for ln = 1:NL,
          if verbose > 0,
            fprintf('\n');
            fprintf('Loop %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',ln,NL,id,File.NT(junc(ln,1)).Number,File.NT(junc(ln,2)).Number,junc(ln,2)+1-junc(ln,1));
          end
          Node(jn).nextnode(ln) =  length(Node)+1;
          Node = pMakeNodes(File,junc(ln,1),junc(ln,2),Interact,Node,n);
%fprintf('List of next nodes:\n');
%for hhh = 1:length(Node),
%  if ~isempty(Node(hhh).nextnode),
%    fprintf('%d ', Node(hhh).nextnode(1));
%  end
%end
%fprintf('\n');
          n = length(Node);
        end
%fprintf('Child nodes of this junction node:\n');
%[jn 0 Node(jn).nextnode]

        return                                  % nothing left to do!

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

        if verbose > 0,
          fprintf('%3d Junction Cluster %4s %4s %4s %4s %4s %4s\n', n, File.NT(r).Number, File.NT(rrr).Number, File.NT(sss).Number, File.NT(ttt).Number, File.NT(uuu).Number, File.NT(u).Number);
        end

        r = rrr;
        s = sss;
        t = ttt;
        u = uuu;

        Node(n).nextnode(1) =  n+1;          % index of next node in tree
        Node = pMakeNodes(File,r,s,Interact,Node,n);

        Node(n).nextnode(2)  = length(Node)+1;
        Node = pMakeNodes(File,t,u,Interact,Node,length(Node));

      end                                  % junction cluster

      Node(n).P    = [0.05*ones(17,1) 0.95*ones(17,1)];
                                            % state to state transitions
      Node(n).PIns	= [0.05 0.95];   % when no previous state

      EndLoop = 1;

    else                                   % not a junction

      % ---------------------------------- Identify basepair or cluster

%[File.NT(a).Number ' ' File.NT(b).Number ' ' File.NT(B).Number]
%[a A AA b B BB]
%G(a,b)
%G(a,B)
%G(a,B-cdepth:B-1)
%G(a+1:a+cdepth,B)

      if ((G(a,b) > 0) && sum(sum(G(a,b-cdepth:b-1))) == 0 ...
                       && sum(sum(G(a+1:a+cdepth,b:B))) == 0), 
                      % a and b interact, but not also with other nearby bases

        % set up basepair
        n=n+1;  
        Node(n).type        = 'Basepair';        % node type
        Node(n).nextnode    = n+1;               % index of next node in tree
        Node(n).LeftLetter  = File.NT(a).Base;
        Node(n).RightLetter = File.NT(b).Base;
        Node(n).Edge        = Interact{a}.Categ(1);
        Node(n).Delete      = 0.01 + n/10000;        % deletion prob
        Node(n).lpar        = [.01*ones(16,1); 0]; % left insertion param
        Node(n).rpar        = [.01*ones(16,1); 0]; % right insertion param
        Node(n).LeftIndex   = a;
        Node(n).RightIndex  = b;
        Score        = pIsoScore(Node(n).Edge,Node(n).LeftLetter, ...
                                 Node(n).RightLetter,Node(n).Delete);
        Node(n).P    = ones(17,1)*Score;         % state to state transitions
        Node(n).PIns = Score;                    % when no previous state
        L = Node(n).lpar(1,1);
        R = Node(n).rpar(1,1);
        X = 0:10;     
        Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
                  * sum(R.^X*exp(-R)./factorial(X));
        Node(n).leftLengthDist  = subPoisson(Node(n).lpar(1,1));
        Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
        Node(n).rightLengthDist = subPoisson(Node(n).rpar(1,1)); 
        Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

        if verbose > 0,
          fprintf('%3d Basepair %4s %4s %s%s %s\n',n, File.NT(a).Number, File.NT(b).Number,File.NT(a).Base,File.NT(b).Base,zEdgeText(File.Edge(a,b)));
        end
        a = a + 1;                                % current base on left
        B = b - 1;

      else 

        % set up cluster node of an appropriate size
        % Later: identify clusters made up entirely of bases on the left
        % or bases on the right.  Now, it combines them, even if they could
        % be left distinct.  Build from the left and from the right, and
        % when they overlap, coalesce them.
%[a b B cdepth]
        b = B;                                 % current base on right

        amax = min(a+cdepth,floor((a+b-0.5)/2));   % how far to look on left
        bmin = max(b-cdepth,ceil((a+b)/2));    % how far to look on right

        X = full(triu(G(a:amax,a:amax)));      % interactions on left
        Y = full(triu(G(bmin:b,bmin:b)));      % interactions on right
        Z = full(G(a:amax,bmin:b));            % interactions across
%  full(X)
%  full(Y)
%  full(Z)
        [s,t] = size(Z);                       % X is s x s, Y is t x t

        ssa = max(find(X(1,:)));               % depth from a on left
        ssb = max(find(Z(:,t)));               % depth from a on right
        if isempty(ssa), ssa = 1; end
        if isempty(ssb), ssb = 1; end
        ss = max(ssa, ssb);

        tta = min(find(Z(1,:)));               % depth from b on left
        ttb = min(find(Y(:,t)));               % depth from b on right
        if isempty(tta), tta = t; end
        if isempty(ttb), ttb = t; end
        tt = min(tta, ttb);
%  [ss tt 2242]
        while sum(sum(Z(1:ss,1:tt-1))) > 0 || ...
              sum(sum(Z(ss+1:s,tt:t))) > 0 || ...
              sum(sum(X(1:ss,ss+1:s))) > 0 || ...
              sum(sum(Y(tt:t,1:tt-1))) > 0,
          ssa = max(find(sum(X(1:ss,:),1)));
          ssb = max(find(sum(Z(:,tt:t),2)));
          if isempty(ssa), ssa = 1; end
          if isempty(ssb), ssb = 1; end
          ss = max(ssa, ssb);

          tta = min(find(sum(Z(1:ss,:),1)));
          ttb = min(find(sum(Y(:,tt:t),2)));
          if isempty(tta), tta = t; end
          if isempty(ttb), ttb = t; end
          tt = min(tta, ttb);
%  [ss tt 3678]
        end

        aa = a - 1 + ss;                      % left extent of cluster
        bb = b - (t - tt);                    % right extent of cluster

        n=n+1;
%[length(Node) n]
%cat(2,Node(:).nextnode)
        Node(n).Delete       = 0.01;
%Node(n)
        Node(n).type         = 'Cluster';        % node type
        Node(n).nextnode     = n+1;              % index of next node in tree
        Node(n).LeftIndex    = [a:aa];           % full range of indices
        Node(n).RightIndex   = [bb:b];

        AllIndices = [a:aa bb:b];

        if verbose > 0,
          fprintf('%3d Cluster  %s:%s %s:%s\n', n, File.NT(a).Number, File.NT(aa).Number, File.NT(bb).Number, File.NT(b).Number);
        end

        % identify which bases are interacting on left and right

        X = X(1:ss,1:ss);                     % reduce to this cluster only
        Y = Y(tt:t,tt:t);
        Z = Z(1:ss,tt:t);

        zs = find(sum(Z,2)');      % bases on left int'ing w/ right
        zt = find(sum(Z,1));       % bases on right int'ing w/ left

        xs = find(sum(X+X',1));    % left with left
        yt = find(sum(Y+Y',1));    % right with right

        leftinter = union(zs,xs);            % list of interacting on left
        leftnum   = [];
        leftnum(leftinter) = 1:length(leftinter); 
                                  % sequential numbering of interacting bases
%leftnum

        rightinter = union(zt,yt);
        rightnum   = [];
        rightnum(rightinter) = 1:length(rightinter);
                                  % sequential numbering of interacting bases
        rightnum(find(rightnum)) = rightnum(find(rightnum))+length(leftinter);
%rightnum

        % list insertion possibilities and corresponding probabilities

%[size(zt) size(yt)]
%zt
%yt
%union(zt,yt)
%Node(n)

        Node(n).Left(1,:)  = union(zs,xs);          % which bases interact
        Node(n).Right(1,:) = union(zt,yt);          % which bases interact
        Node(n).PIns = [0.00001 0.99999];           % when no previous state
        Node(n).P    = ones(17,1) * Node(n).PIns;

        % create a list of insertions according to what is observed here

        e = [Node(n).Left(1,:) Node(n).Left(1,end)+Node(n).Right(1,:)];
                                                    % bases used, left to right

        d = diff(e);                                % diffs in positions used
        h = find(d>1);                              % where insertions occur
        cc = 1;
        for aaa = 1:length(h),
          Node(n).Insertion(cc).Position = h(aaa);
          Node(n).Insertion(cc).LengthDist = subPoisson(d(h(aaa))-1);
          Node(n).Insertion(cc).LetterDist = [1 1 1 1]/4;  % WRONG!!!

          if verbose > 0,
            fprintf('    %d insertions between %s%s and %s%s\n', ...
            d(h(aaa))-1, File.NT(AllIndices(h(aaa))).Base, ...
                         File.NT(AllIndices(h(aaa))).Number, ...
                         File.NT(AllIndices(h(aaa)+d(h(aaa)))).Base, ...
                         File.NT(AllIndices(h(aaa)+d(h(aaa)))).Number);
          end
          cc = cc + 1;
        end

        % create a list of interacting bases
        % interactions between left and right
        K = 1;                                % counter
        [i,j] = find(Z);                      % interacting pairs
        for k = 1:length(i),                  % loop through them
            Node(n).IBases(K,:) = [leftnum(i(k)) rightnum(j(k))];
            i1 = i(k) + a - 1;             % index of first base
            i2 = j(k) + bb - 1;             % index of second

            if verbose > 0,
              fprintf('     %s %s %s %5.1f\n', File.NT(i1).Number, zEdgeText(File.Edge(i1,i2)), File.NT(i2).Number, full(File.Edge(i1,i2)));
%              fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
            end

            Node(n).Score(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
                                   File.NT(i1).Code, File.NT(i2).Code);

            Node(n).SubsProb(K,:) = pIsoScore(File.Edge(i1,i2), ...
                                   File.NT(i1).Code, File.NT(i2).Code,0.0);

            K  = K + 1;
        end

        % interactions between left and left
        [i,j] = find(X);                     % interacting pairs
        for k = 1:length(i),                 % loop through them
            Node(n).IBases(K,:) = [leftnum(i(k)) leftnum(j(k))];
            i1 = i(k) + a - 1;             % index of first base
            i2 = j(k) + a - 1;             % index of second

            if verbose > 0,
              fprintf('     %s %s %s %5.1f\n', File.NT(i1).Number, zEdgeText(File.Edge(i1,i2)), File.NT(i2).Number, File.Edge(i1,i2));
%fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
            end

            Node(n).Score(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
                                   File.NT(i1).Code, File.NT(i2).Code);
            K  = K + 1;
        end

        % interactions between right and right
        [i,j] = find(Y);                      % interacting pairs
        for k = 1:length(i),                  % loop through them
            Node(n).IBases(K,:) = [rightnum(i(k)) rightnum(j(k))];
            i1 = i(k) + bb - 1;             % index of first base
            i2 = j(k) + bb - 1;             % index of second

            if verbose > 0,
              fprintf('     %s %s %s %5.1f\n', File.NT(i1).Number, zEdgeText(File.Edge(i1,i2)), File.NT(i2).Number, File.Edge(i1,i2));
%fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
            end

            Node(n).Score(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
                                   File.NT(i1).Code, File.NT(i2).Code);
            K  = K + 1;
        end

%Node(n).IBases

%pause

        a = aa + 1;                           % current base on left
        B = bb - 1;                           % current base on right

      end                                          % basepair or cluster

      % ------------------- check for hairpin and if not, probe for insertions

      if (a == B) || (sum(sum(G(a:B,a:B))) == 0),
        n = n + 1;
        Node(n).type = 'Hairpin';
        Node(n).MiddleIndex = [a:B];
        if a > B,
          Node(n).MiddleIndex = [a];
        end
        Node(n).P       = ones(17,1);
        Node(n).PIns    = 1;
        Node(n).subtype = 'XXXX';                  % revise this later!
        if verbose > 0,
          fprintf('%3d Hairpin  %s:%s\n', n, File.NT(a).Number, File.NT(B).Number);
        end
        EndLoop = 1;

      else                                 % probe for insertions

        LeftIns = 0;
 
        aa = a;                            % initial values of these
        BB = B;

        while sum(G(a,(a+1):B)) == 0,      % no interaction w/in this loop
          if verbose > 0,
            fprintf('    Insertion %3s\n', File.NT(a).Number);
          end
          LeftIns = LeftIns + 1;           % increase mean number of insertions
          a = a + 1;                               % next base on left
        end

        RightIns = 0;

        while sum(G(a:(B-1),B)) == 0,      % no interaction w/in this loop
          if verbose > 0,
            fprintf('    Insertion     %4s\n', File.NT(B).Number);
          end
          RightIns = RightIns + 1;         % increase mean number of insertions
          B = B - 1;                               % next base on right
        end
 
        if LeftIns > 0 || RightIns > 0,
          if strcmp(Node(n).type,'Basepair'),   % add insertions to basepair
            Node(n).lpar = LeftIns  * [ones(16,1); 0];
            Node(n).rpar = RightIns * [ones(16,1); 0];
            Node(n).leftLengthDist  = subPoisson(LeftIns);
            Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
            Node(n).rightLengthDist = subPoisson(RightIns); 
            Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];
          elseif strcmp(Node(n).type,'Cluster'),
            n = n + 1;                          % initial node after cluster
            Node(n).type = 'Initial';
            Node(n).nextnode  = n+1;            % index of next node in tree
            Node(n).lpar = LeftIns;
            Node(n).rpar = RightIns;
            Node(n).LeftIndex  = aa;
            Node(n).RightIndex = BB;
            L = Node(n).lpar(1,1);
            R = Node(n).rpar(1,1);
            X = 0:10;
            Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
                      * sum(R.^X*exp(-R)./factorial(X));%normalization constant
            Node(n).leftLengthDist  = subPoisson(Node(n).lpar);
            Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
            Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
            Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

            if verbose > 0,
              fprintf('%3d Initial %s:%s and %s:%s\n', n, File.NT(aa).Number, File.NT(aa+LeftIns).Number, File.NT(BB).Number, File.NT(BB+RightIns).Number);
            end
          end
        end
      end                                 % hairpin or insertions
    end                                    % junction and junction cluster

  else
    fprintf('Found a case with a > b  ========================\n')
    break
  end                                   % if (a < b)

end                                       % while loop

% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
