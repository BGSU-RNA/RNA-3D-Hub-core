% pMakeNodes(File,NTNumber) makes a secondary structure node model based on
% the Edge interaction matrix in File, starting at NTNumber and ending at
% LastNTNumber

function [Node] = pMakeNodes(File,NTNumber,LastNTNumber,Interact,Node,n)

cdepth  = 6;                      % how far to look ahead for a cluster
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
G = (E < 15) .* (E ~= 0);                  % consider basepairing only

if nargin == 2,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
end

if nargin < 4,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
  [y,L] = sort(E(a,k));                    % sort by edge interaction category
  Interact{a}.Categ = File.Edge(a,k(L));   % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

a  = Indices(1);                           % first index; current index
A  = a;                                    % previous interacting base on left
B  = LastNTNumber;                         % next base on right
AA = a;                                    % previous cWW base on left
BB = a;                                    % previous cWW base on right

fprintf('Loop %4s %4s\n', File.NT(a).Number, File.NT(B).Number);

% Initial node creation ------

n = n+1;                                   % move to next node
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).lpar      = 0;                     % left insertion parameter
Node(n).rpar      = 0; 
Node(n).Bl        = 0;                     % actual left base
Node(n).Br        = 0;
                                           % probe for insertions
while length(Interact{a}.Categ) == 0,      % no interaction
  Node(n).lpar = Node(n).lpar + 1;         % increase mean number of insertions
  a = a + 1;                               % move to next base
end

% Later: deal with the range from b to N

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

      if (sum(sum(G(t:u,t:u) == 1))         > 0) & ...
         (sum(sum(G(a+1:b-1,a+1:b-1) == 1)) > 0),
           % there are helices between a and b and between b+1 and B
           n = n+1;                        % move to next node

           C =   sum(sum(G(r:r+jcdepth,t:t+jcdepth))) ... % junction cluster inter
               + sum(sum(G(r:r+jcdepth:u-jcdepth:u))) ...
               + sum(sum(G(s-jcdepth:s,t:t+jcdepth))) ...
               + sum(sum(G(s-jcdepth:s,u-jcdepth:u)));

           if C == 0,
             Node(n).type        = 'Junction';       % plain junction
           else
             Node(n).type        = 'JunctionCluster';  % 

             % find extent of interactions between loops
             % first, extent of interactions between loops

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

cat(2,File.NT([r rrr sss s t ttt uuu u]).Number)

             Node(n).LeftIndices   = [r:rrr];
             Node(n).MiddleIndices = [sss:ttt];
             Node(n).RightIndices  = [uuu:u];

             r = rrr;
             s = sss;
             t = ttt;
             u = uuu;

           end                                  % junction cluster

           Node(n).nextnode(1) =  n+1;          % index of next node in tree
           Node = pMakeNodes(File,r+1,s-1,Interact,Node,n);
           Node(n).nextnode(2)  = length(Node)+1;
           Node = pMakeNodes(File,t+1,u-1,Interact,Node,length(Node));

           Node(n).P    = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                 % state to state transitions
           Node(n).PIns	= [0.05 0.95];   % when no previous state

           EndLoop = 1;
      end

      % ---------------------------------- Identify basepair or cluster

      if ((G(a,B) > 0) && sum(sum(G(a,B-cdepth:B-1))) == 0 ...
                         && sum(sum(G(a+1:a+cdepth,B))) == 0), 
                    % a and B interact, but not also with other nearby bases

          % set up basepair
          n=n+1;  
          Node(n).type        = 'Basepair';        % node type
          Node(n).nextnode    = n+1;               % index of next node in tree
          Node(n).LeftLetter  = File.NT(a).Base;
          Node(n).RightLetter = File.NT(b).Base;
          Node(n).Edge        = Interact{a}.Categ(1);
          Node(n).Delete      = 0.1;                  % deletion prob
          Node(n).lpar        = [.001*ones(16,1); 0]; % left insertion param
          Node(n).rpar        = [.001*ones(16,1); 0]; % right insertion param
          Node(n).LeftIndex   = a;
          Node(n).RightIndex  = b;
          Node(n).Bl          = [str2num(File.NT(a).Number)]; %nucleotide number
          Node(n).Br          = [str2num(File.NT(b).Number)]; % nucl number 

          a = a + 1;                                % current base on left
          B = B - 1;                                % current base on right

      else 

             % set up cluster node of an appropriate size

             b = B;                                 % current base on right

             amax = min(a+cdepth,b-1);              % how far to look on left
             bmin = max(amax+1,b-cdepth);               

 fprintf('Cluster %4s %4s %4s %4s\n', File.NT(a).Number, File.NT(amax).Number, File.NT(bmin).Number, File.NT(b).Number);

             X = triu(G(a:amax,a:amax));            % interactions on left
             Y = triu(G(bmin:b,bmin:b));            % interactions on right
             Z = G(a:amax,bmin:b);                  % interactions across
full(X)
full(Y)
full(Z)
             [s,t] = size(Z);                       % X is s x s, Y is t x t
             ss = max( max(find(Z(:,t))), max(find(X(1,:)))); % depth from a
             if isempty(ss), ss = 1; end
             tt = min( min(find(Z(1,:))), min(find(Y(:,t)))); % depth from b
             if isempty(tt), tt = t; end
[ss tt 0]
             while sum(sum(Z(1:ss,1:tt-1))) > 0 || ...
                   sum(sum(Z(ss+1:s,tt:t))) > 0 || ...
                   sum(sum(X(1:ss,ss+1:s))) > 0 || ...
                   sum(sum(Y(tt:t,1:tt-1))) > 0,
               ss = max(max(find(sum(Z(:,tt:t)')')),max(find(sum(X(1:ss,:)))));
               tt = min(min(find(sum(Z(1:ss,:)))),min(find(sum(Y(:,tt:t)')')));
[ss tt 0]
             end

             aa = a - 1 + ss;                      % left extent of cluster
             bb = b - (t - tt);                 % right extent of cluster

             Node(n).LeftIndices   = [a:aa];
             Node(n).RightIndices  = [bb:b];

             a = aa + 1;                           % current base on left
             B = bb;                               % last base accounted for

      end                                    % basepair or cluster

      % ------------------- check for hairpin and if not, probe for insertions

 fprintf('%3d %s %4s %4s\n',n, Node(n).type, File.NT(a-1).Number, File.NT(B+1).Number);

size(a)
size(B)
[a B]

      if (a == B) || (sum(sum(G(a:B,a:B))) == 0),
        n = n + 1;
        Node(n).type = 'Hairpin';
        Node(n).Indices = [aa:B];
        EndLoop = 1;
      else
          LeftIns = 0;
 
          while length(Interact{a}.Categ) == 0,      % no interaction
            LeftIns = LeftIns + 1;    % increase mean number of insertions
            a = a + 1;                               % next base on left
          end

          RightIns = 0;

          while length(Interact{B}.Categ) == 0,      % no interaction
            RightIns = RightIns + 1;  % increase mean number of insertions
            B = B - 1;                               % next base on right
          end

          if strcmp(Node(n).type,'Basepair'),
            Node(n).lpar = LeftIns  * [ones(16,1); 0];
            Node(n).rpar = RightIns * [ones(16,1); 0];
          elseif strcmp(Node(n).type,'Cluster'),
            n = n + 1;
            Node(n).type = 'Initial';
            Node(n).lpar = LeftIns;
            Node(n).rpar = RightIns;
          end

      end                                 % hairpin or insertions
    else
      fprintf('Found a case with a > b  ========================\n')
    end                                   % if (a < b)
end                                       % while loop

