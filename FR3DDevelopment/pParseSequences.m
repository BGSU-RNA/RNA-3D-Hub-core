% pParseSequences(node,Sequence,NumSequences) uses the model in Node to
% parse sequences First to Last in Sequence

% The program makes plans to use previous computations to speed up parsing
% of later sequences, so it is important to start with sequence 1 and proceed
% one after another through Sequence.

% Reuse = 1 means save data for use with later sequences

function [Node,Sequence] = pParseSequences(Node,Sequence,First,Last,Reuse)

if nargin < 5,
  Reuse = 1;
end

%clf

K = length(Sequence);
N = length(Node);           % number of nodes, including initial and hairpin

% ----------------------------------------- Internal calculations

% These should be separated out into a new program, since they may change
% the length of Sequence

if First == 1,                            % before parsing any sequences ...

  Sequence(1).EstTime = 0;                % create this field

  i = [];

  % --------------- exclude sequences with letters other than ACGU.-

  for k=1:length(Sequence),
    if min(ismember(Sequence(k).X,'ACGU.-')) == 1,  % only these letters
      i = [i k];
    end
  end
 
  Sequence = Sequence(i);                       % remove those w/ difft letters
  K = length(Sequence);

  for k = 1:K,
    Seq{k} = Sequence(k).X;                     % for below
  end

  if Reuse == 1,

    % ---------------------------------- look for overlapping subsequences

    Over = sparse(zeros(K));

    DeleteAfter = 1:K;                            % when to delete parse info
    for k = 2:K,
      Sequence(k).Overlap = pLCS(Sequence(k).X, Seq(1:(k-1)),110);
      for r = 1:length(Sequence(k).Overlap),      
                                   % for each earlier overlapping sequence,
        DeleteAfter(Sequence(k).Overlap(r).Seqnum) = k; 
                                   % wait until seq k to delete
        Over(k,Sequence(k).Overlap(r).Seqnum) = 1;
      end
    end

    for k=1:K,
      Sequence(k).DeleteAfter = DeleteAfter(k);
    end

    fprintf('Analyzed sequences for redundant subsequences\n');
  else

    for k=1:K,
      Sequence(k).Overlap     = [];
      Sequence(k).DeleteAfter = k;
    end

  end

  % ---------------------------------- Precompute logarithms, etc. for Node

  pairs = {'AA', 'CA', 'GA', 'UA', 'AC', 'CC', 'GC', 'UC', 'AG', 'CG', 'GG', 'UG', 'AU', 'CU', 'GU', 'UU', 'DEL'};

  k = 10;                     % maximum number of insertions to calculate for

  warning off

  for n=1:length(Node),
    switch Node(n).type,
      case {'Initial','Junction','JunctionCluster','Hairpin','Alternative'}
        Node(n).numstates = 1;
      case {'Basepair','Cluster'}
        Node(n).numstates = 2;    % deleted or active
    end

    switch Node(n).type,
      case {'Basepair','Hairpin','Junction','JunctionCluster','Cluster','Alternative'}
        [s,t] = size(Node(n).P);
        Node(n).lP = log(Node(n).P);
        if t == 17,                                 % add one more for .A
          Node(n).lP = [Node(n).lP -Inf*ones(s,1)];
        end
        Node(n).lPIns      = [log(Node(n).PIns) -Inf];           
      case 'Initial'
        Node(n).lP = zeros(17,1);                   % for transition to Initial
        Node(n).lPIns = zeros(17,1);
    end
    
    switch Node(n).type,
      case {'JunctionCluster'}
        Node(n).LLIP       = log(Node(n).LIP);  
        Node(n).LMIP       = log(Node(n).MIP);
        Node(n).LRIP       = log(Node(n).RIP);           
    end
    
    switch Node(n).type,
      case {'Cluster'}
        Node(n).LLIP       = log(Node(n).LIP);           
        Node(n).LRIP       = log(Node(n).RIP);           
    end

    switch Node(n).type,
      case {'Initial','Basepair'}
        for s=1:length(Node(n).lpar(:,1)),
          l = Node(n).lpar(s);
          Node(n).lip(s,:) = exp(-l)*((ones(1,k+1)*l).^(0:k))./[1 cumprod(1:k)];
          r = Node(n).rpar(s);
          Node(n).rip(s,:) = exp(-r)*((ones(1,k+1)*r).^(0:k))./[1 cumprod(1:k)];
        end
        Node(n).loglip = log(Node(n).lip)-log(Node(n).Z);  % the Z(n) is the normalizing constant
        Node(n).logrip = log(Node(n).rip);
    end
  end

%   warning on

  fprintf('Did initial internal calculations\n');

end

% ------------------------------ Estimate parsing time

if (First >= 2) & (Sequence(1).EstTime == 0),    % after first is parsed
  Sequence(1).EstTime = Sequence(1).time;

  for k=2:K,
    Done = zeros(length(Sequence(k).X));
    for r=1:length(Sequence(k).Overlap),
      a = Sequence(k).Overlap(r).XStart:Sequence(k).Overlap(r).XStop;
      Done(a,a) = ones(length(a),length(a));
    end
    Sequence(k).EstTime = Sequence(1).time*sum(sum(1-Done))/length(Sequence(1).X)^2;
  end
else
  Sequence(1).EstTime = 0;
end

for k=First:K,
  Sequence(k).time = 0;
end

% -------------------------------------------- Loop through sequences

for k=First:min(Last,K),
  Sequence(k).time = cputime;
  L = length(Sequence(k).X);

  % ------------------------------------- initialize parsing information

  for n=1:N,
    for s=1:Node(n).numstates,
      Sequence(k).P(n,s).mp   = -Inf*ones(L,L);% maximum probability matrix
      Sequence(k).P(n,s).sub  = uint16(zeros(L,L)); % indices explained by next node
      Sequence(k).P(n,s).sub2 = uint16(zeros(L,L)); % indices explained by next node
      Sequence(k).P(n,s).rmi  = uint16(zeros(L,L));
      Sequence(k).P(n,s).next = uint16(zeros(L,L)); % state and number of next node
      Sequence(k).Computed    = uint16(zeros(L,L));
      Sequence(k).TraceInfo   = [];
    end
  end

  % ------------------------------------- fill in data from previous sequences

  if k > 1,
    for r=1:length(Sequence(k).Overlap),
      a = Sequence(k).Overlap(r).XStart:Sequence(k).Overlap(r).XStop;
      b = Sequence(k).Overlap(r).YStart:Sequence(k).Overlap(r).YStop;
      c = Sequence(k).Overlap(r).Seqnum;
      
      Sequence(k).Computed(a,a)    = ones(length(a),length(a));

      for n=1:length(Node),
        for s=1:Node(n).numstates,
          Sequence(k).P(n,s).mp(a,a)   = Sequence(c).P(n,s).mp(b,b);
          Sequence(k).P(n,s).sub(a,a)  = uint16(double(Sequence(c).P(n,s).sub(b,b))+double(a(1))-double(b(1)));
          Sequence(k).P(n,s).next(a,a) = Sequence(c).P(n,s).next(b,b);
          Sequence(k).P(n,s).sub2(a,a) = uint16(double(Sequence(c).P(n,s).sub2(b,b))+double(a(1))-double(b(1)));
          Sequence(k).P(n,s).rmi(a,a)  = uint16(double(Sequence(c).P(n,s).rmi(b,b))+double(a(1))-double(b(1)));
        end
      end
    end
  end

  % -------------------------------------- Parse this sequence

  Sequence(k) = pParseOneSequence(Node,Sequence(k));

  Sequence(k).time = (cputime-Sequence(k).time)/60;


  % -------------------------------------- Delete unneeded computations

  for j=1:k,
    if Sequence(j).DeleteAfter <= k,
%      Sequence(j).P = [];
%      Sequence(j).Computed = [];
    end
  end

  % -------------------------------------- Display computation time
  
  Total = sum(cat(1,Sequence(1:k).time));
  s = whos('Sequence');                     % size of stored sequence data

  fprintf('Sequence %s took %5.3f min. Total %5.3f min, %5.3f per sequence. Memory %6.2f MB\n', Sequence(k).FastaNum, Sequence(k).time, Total, Total/k, s.bytes/(1024^2));

  drawnow

end

% [(1:K)' cat(1,Sequence(:).EstTime) cat(1,Sequence(:).time) cat(1,Sequence(:).DeleteAfter) ] % display parse times




