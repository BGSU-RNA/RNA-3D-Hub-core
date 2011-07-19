% pMakeNodes(File,NTNumber,LastNTNumber,Truncate,Interact,Node,n) makes a secondary structure node model based on the Edge interaction matrix in File, starting at NTNumber and ending at LastNTNumber.  It assigns various nodes consistent with this secondary structure.  Truncate indicates where to put * hairpins.  Interact, Node, and n are optional parameters specified when pMakeNodes is called by itself.

function [Node] = pMakeNodes(File,Param,NTNumber,LastNTNumber,Truncate,Data,Node,n)

if nargin < 2,
  Verbose = 1;
end

if nargin < 3,
  NTNumber = 1;
end

Verbose         = Param(1);
method          = 4;             % method for assigning pair subst probs
Extension       = 1;             % whether to extend stems with no LR inter
AdjustSubsForLR = 1;             % adjust ins, basepair subs probs for LR inter
cdepth          = 10;            % how far to look ahead for a cluster
usenear         = 0;             % use near pairs for basepair probabilities

jcdepth         = 4;             % how far to look for a junction cluster

% Parameters stored in Param:
% Param(1) verbose
% Param(2) method to use for basepair isostericity
% Param(3) recognize extensible helices and model them as such
% Param(4) adjust substitution probabilities for long-range interactions
% Param(5) how far to look ahead for local basepair interactions
% Param(6) use near interactions


if length(Param) > 5,
  usenear = Param(6);
end

if length(Param) > 4,
  cdepth = Param(5);
end

if length(Param) > 3,
  AdjustSubsForLR = Param(4);
end

if length(Param) > 2,
  Extension = Param(3);
end

if length(Param) > 1,
  method  = Param(2);
end

if nargin < 5,
  Truncate = [];
end

if nargin < 8,
  n=0;                            % current node number
end

% -------------------------- if File is a text string (filename), load the file

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% ----------------- if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char'),
  NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell'),
  NTNumber = zIndexLookup(File,NTNumber);
end

% ----------------- if Truncate is a cell array of numbers, look up the indices

if strcmp(class(Truncate),'char'),
  Truncate = {Truncate};
end

if strcmp(class(Truncate),'cell'),
  if isempty(Truncate),
    Truncate = [];
  else
    Truncate = zIndexLookup(File,Truncate);
  end
end

% ------------------------------------------ Set key variables

N = length(File.NT);                       % number of nucleotides in File
DelProb = 0.01;                            % nominal deletion probability
                                           % for basepairs
TertiaryFreeNode = 0;                      % first node in this stem making
                                           % no tertiary intearctions beyond it

if ~isfield(File,'BasePhosphate'),
  File.BasePhosphate = sparse(zeros(N,N));
end

if ~isfield(File,'BaseRibose'),
  File.BaseRibose = sparse(zeros(N,N));
end

if nargin < 4,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
elseif strcmp(class(LastNTNumber),'char'),
  LastNTNumber = zIndexLookup(File,{LastNTNumber});
end



%NTNumber
%LastNTNumber
%full(File.Edge(NTNumber:(NTNumber+10),(LastNTNumber-10):LastNTNumber))



% ------------------------------------------ Store indices of interacting bases

  load PairExemplars

if nargin < 6,

  E = abs(fix(File.Edge));                   % don't distinguish subcategories
  G = E .* (E < 16) .* (E ~= 0);             % consider basepairing only,
                                             % including bifurcated, water-ins

  if usenear > 0,
    G = G +  E .* (E > 100);                 % near basepairs too
  end

  H = (G ~= 0) .* max(File.Crossing == 0, abs(G) == 1) ;
                                             % 1 for nested pairs, 0 otherwise

  J = abs(G .* (File.Crossing >  0));        % long-range basepairs only

  GG = G .* (abs(fix(G)) ~= 9);              % eliminate cSH pairs for hairpins
  GG = G;
  for i = 1:(length(GG(:,1))-1),
    GG(i,i+1) = 0;                           % eliminate pairs btw adjacent
    GG(i+1,i) = 0;
  end

  for a = 1:N,                             % loop through nucleotides
   k = find(G(a,:));                       % find indices of interacting bases
   [y,L] = sort(E(a,k));                   % sort by edge interaction category
   Interact{a}.Categ = abs(File.Edge(a,k(L)));   % store categories
   Interact{a}.Index = k(L);               % store indices of interacting bases
  end

  % ------------------------------------------ prepare to identify motifs

  HasMotif     = zeros(1,length(File.NT));
  HasGUPacking = zeros(1,length(File.NT));   % highly-conserved motif

  if isfield(File,'Nucl'),
    for i = 1:length(File.NT),
      if ~isempty(File.Nucl(i).Motif),
        HasMotif(i) = 1;
        for m = 1:length(File.Nucl(i).Motif),
          if ~isempty(strfind(File.Nucl(i).Motif(m).Name,'GU_packing')),
            HasGUPacking(i) = 1;
          end
        end
      end
    end
  end

  Data.Interact = Interact;
  Data.HasMotif = HasMotif;
  Data.HasGUPacking = HasGUPacking;
  Data.E = E;
  Data.G = G;
  Data.H = H;
  Data.J = J;
else
  Interact = Data.Interact;
  HasMotif = Data.HasMotif;
  HasGUPacking = Data.HasGUPacking;
  E            = Data.E;
  G            = Data.G;
  H            = Data.H;
  J            = Data.J;
end

% ------------------------------------------ Set up initial values of counters
a  = NTNumber;                             % first index; current index
A  = a;                                    % previous interacting base on left
AA = a;                                    % previous cWW base on left

B  = LastNTNumber;                         % next base on right
BB = a;                                    % previous cWW base on right


if Verbose > 0,
  fprintf('Loop %4s %4s\n', File.NT(a).Number, File.NT(B).Number);
end

% Initial node creation -------------------------------------------------

n = n+1;                                   % move to next node

pMakeNodesNewNode;                         % set up blank node with all fields
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).LeftIndex = a;                     % index of first base on left
Node(n).RightIndex= B;                     % index of first base on right


pMakeNodesProbeForInsertions;              % probe for insertions, each strand

% ---------------------------------------------------------------------------

EndLoop = 0;                               % flag for the end of the loop

while (EndLoop == 0) && (a <= LastNTNumber), % while not the end of the loop,

    % ---------------------------------- Check for junction

    % check to see if a is now in a new nested loop

    r = a;                               % leftmost index of a cWW
    rr= a;                               % start of current loop

    while sum(H(r,(r+1):B)) == 0 && r < B, % if r does not make a nested pair,
      r = r + 1;
    end
    
    s = Interact{r}.Index(1);            % what it interacts with
    t = s+1;                             % next after that
    u = B;                               % end of current known loop

    if (sum(sum(G(t:u,t:u)==1)) > 0) && (sum(sum(G(r:s,r:s)==1)) > 0),
            % there are nested cWW pairs between r and s and between t and u
            % use cWW pairs to avoid pairs between i and i+1 making junctions

      if Verbose > 1,
        fprintf('Found nested interactions between %s and %s and between %s and %s\n', File.NT(r).Number, File.NT(s).Number, File.NT(t).Number, File.NT(u).Number);
      end

      pMakeNodesJunction                   % make models for junctions
      return                               % nothing left to do!

    else                                   % not a junction

      % ---------------------------------- Identify basepair or cluster

      aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
      BBB = max([1 B-cdepth ceil((a+B)/2)]);
      LS = (a+1):aaa;                         % left strand
      RS = BBB:(B-1);                         % right strand

      if HasMotif(a),   % ---------------------- Insert motif model when needed
 
        pMakeNodesMotif

      elseif (H(a,B) > 0) && ((cdepth == 0) || ...
        (sum(sum(G(a,[LS RS]))) == 0 && sum(sum(G([LS RS],B))) == 0)), 
                      % a and B interact, but not also with other nearby bases
        pMakeNodesBasepair                       % add basepair with insertions

      else     % a and B also interact with nearby bases - use a cluster node

        % [a aaa BBB B]
        % zShowInteractionTable(File,unique([a:aaa BBB:B]));

        if abs(B - a) > 7,                         % far enough from hairpin
          pMakeNodesCluster;
        else
          pMakeNodesHairpin;
          EndLoop = 1;
          if Verbose > 0,
            disp('Making a hairpin instead of inserting a cluster');
          end
        end
        % pause

      end                                          % basepair or cluster

      % ------------------- check for tertiary interactions in this stem

      pMakeNodesCheckForExtensibility

      % ------------------- check for truncation and hairpin

      if (EndLoop == 0),
        if ismember(a,Truncate) || ismember(a-1,Truncate) || isempty(File.NT(a).Base),

          pMakeNodesTruncate                % add * hairpin

        elseif (a == B) || ((sum(sum(abs(GG(a:B,a:B)))) == 0)), 
                                            % no nucleotides left, or
                                            % no further basepairs except cSH
          if (TertiaryFreeNode > 0) && isempty(Truncate) && Extension > 0,
                     % extensible region, not a truncated model
            pMakeNodesExtraBasepairs;       % add extra basepairs if extensible
          end

          pMakeNodesHairpin
        else

          pMakeNodesProbeForInsertions       % add insertions if needed

        end                                  % hairpin or insertions
      elseif TertiaryFreeNode > 0 && isempty(Truncate) && Extension > 0,
        pMakeNodesExtraBasepairs;       % add extra basepairs if extensible
      end                                    % if EndLoop == 0
    end                                      % junction and junction cluster
end                                       % while (EndLoop == 0) & (a <= N),

% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
