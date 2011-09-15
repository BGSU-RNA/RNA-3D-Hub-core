% pMakeModelFromSearchSaveFile(Search) creates an SCFG/MRF Node variable corresponding to the model in Search

% pMakeModelFromSearchSaveFile('LIB00002 IL 2008-03-20_23_29_25-Sarcin_13_flanked_by_cWW_in_1s72')
% Search = 'LIB00002 IL 2008-03-20_23_29_25-Sarcin_13_flanked_by_cWW_in_1s72';

% load LIB00014_IL_tSH-tSH-tHS-tHS.mat
% pMakeModelFromSearchSaveFile(Search,'IL',1);

function [Node,Truncate] = pMakeModelFromSearchSaveFile(Search,Param)

if nargin < 2,
  Param   = 0;
  Verbose = 0;
else
  Verbose = Param(1);
end

% ----------------------------------- Load Search from filename, if applicable

if strcmp(class(Search),'char'),
  load(['MotifLibrary' filesep Search],'Search','-mat');
end

% ----------------------------------- Gather basic information about the search

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

f = Search.Candidates(:,N+1);           % file numbers of motifs

% the following line assumes that the GU packing motif has the text
% 'GU_packing' in its name and that the GU pair are nucleotides 1 and 2

%if ~isempty(strfind(Search.SaveName,'GU_packing')), % GU packing motif
%  N = 2;                                          % skip last nucleotide
%end

File = Search.File(f(1));                      % file of query motif
NTNumber = double(Search.Candidates(1,1));     % index of first NT
LastNTNumber = double(Search.Candidates(1,N)); % index of last NT

% --------------------------------------- Find locations of truncations

if isfield(Search.Query,'MaxDiffMat'),
  MaxDiff = diag(Search.Query.MaxDiffMat,1);
else
  MaxDiff = Inf*ones(1,N-1);
end

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Search.Candidates(c,1:N))))-1);
end

Truncate = [];

for n = 1:(N-1),
  if (MaxDiff(n) == Inf) || (maxinsert(n) > 5),   % if only few insertions
    Truncate = [Truncate n+1];                    % truncate the strand here
  end
end

Truncate

F.Edge = pConsensusInteractions(Search);


% -------------------------------- Make the model for the consensus structure

i = Search.Candidates(1,1:N);                   % indices of query motif

if Verbose > 0,
  fprintf('Query motif interactions:\n');
  full(fix(File.Edge(i,i)))

  fprintf('Consensus interactions:\n')
  full(F.Edge)                                  % consensus interactions
end

F.NT = File.NT(Search.Candidates(1,1:N));   % use the first candidate as model
F.Crossing = zeros(N,N);                    % small enough, pretend none

% -------------- number strands starting at 1, at 101, at 201, etc.

if length(Truncate) > 0,                    % at least two strands
  b = 1:N;
  for t = 1:N,
    b(t) = b(t) + 100*sum(t >= Truncate);
  end
  binv = 1:max(b);                                  % invert the spreading

size(binv)
size(b)
N

  binv(b) = 1:N;
else
  b = 1:N;
  binv = 1:N;
end

FF.Filename      = File.Filename;
FF.Edge(b,b)     = F.Edge;                        % spread the strands out
FF.NT(b)         = F.NT;
FF.Crossing(b,b) = F.Crossing;

Param(1) = 1;

Node = pMakeNodes(FF,Param,1,b(N),Truncate);          % make the SCFG/MRF model

for n = 1:length(Node),
  Node(n).LeftIndex    = binv(Node(n).LeftIndex);
  Node(n).RightIndex   = binv(Node(n).RightIndex);
  Node(n).MiddleIndex  = binv(Node(n).MiddleIndex);
  Node(n).InterIndices = binv(Node(n).InterIndices);
end

% ---------------------------- Set parameters for the nodes from instances

for n = 1:length(Node),
  switch Node(n).type
  case 'Initial'

    % we should probably allow for stray bases at the beginning

  case 'Basepair'
    a = Node(n).LeftIndex;                   % which NT of the query motif
    b = Node(n).RightIndex;                  % right NT of the query motif


disp('Getting consensus for a basepair');

    Score = pConsensusPairSubstitution(a,b,f,Search.File,F,L,Search,Verbose);

    if Verbose > 0,
      fprintf('Original substitution probabilities\n');
      Node(n).SubsProb

      fprintf('Consensus substitution probabilities\n');
      Score
    end

    Node(n).SubsProb = Score;

    % ----------------------------- tally insertions on the left
    inscount = [];
    letter = [0 0 0 0];                      % record which bases occur
    if n < length(Node),
      for c = 1:L,
        inscount(c) = Node(n+1).LeftIndex(1)-Node(n).LeftIndex(1)-1;
      end
    end

    lld = ones(1,max(inscount)+2);           % Dirichlet distribution
    for c = 1:L,
      lld(inscount(c)+1) = lld(inscount(c)+1) + 1;
    end

    Node(n).leftLengthDist = lld / sum(lld);    % normalize

    % ----------------------------- tally insertions on the right
    inscount = [];
    if n < length(Node),
      for c = 1:L,
        inscount(c) = Node(n).RightIndex(1)-Node(n+1).RightIndex(1)-1;
      end
    end

    rld = ones(1,max(inscount)+2);           % Dirichlet distribution
    for c = 1:L,
      rld(inscount(c)+1) = rld(inscount(c)+1) + 1;
    end

    Node(n).rightLengthDist = rld / sum(rld);    % normalize

    % the following line assumes that the GU packing motif has the text
    % 'GU_packing' in its name and that the GU pair are nucleotides 1 and 2

%    if ~isempty(strfind(Search.SaveName,'GU_packing')),  % GU packing motif
%      P = Node(n).SubsProb;
%      P(3,4) = 0;
%      P = 0.2 * P / sum(sum(P));
%      P(3,4) = 0.8;                             % 80% probability of being GU
%      Node(n).SubsProb = P;
%    end

  case 'Cluster'
    Indices = [Node(n).LeftIndex(Node(n).Left) ...
               Node(n).RightIndex(Node(n).Right)];
    for ii = 1:length(Node(n).IBases(:,1)),
      a = Indices(Node(n).IBases(ii,1));
      b = Indices(Node(n).IBases(ii,2));
      Score = pConsensusPairSubstitution(a,b,f,Search.File,F,L,Search,Verbose);
      Node(n).SubsProb(:,:,ii) = Score;
      if Verbose > 0,
        fprintf('\n');
      end
    end  

  case 'Junction'

  end

  if Verbose > 0,
    fprintf('\n')
  end
end
