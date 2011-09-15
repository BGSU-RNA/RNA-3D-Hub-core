% pMakeNodesExtraBasepairs adds cWW basepairs with relatively high deletion probability, which is used in stems that make no tertiary interactions

% the number of possible extra basepairs is hard to pinpoint!

  nn = n;                               % store current node number
  while fix(abs(Node(n).Edge)) ~= 1,    % work backward through non-cWW pairs
    n = n - 1;
  end

  ClosingPairNode = Node((n+1):nn);     % store for later

  if Verbose > 0,
    fprintf('Found %d non-cWW closing basepair\n', length(ClosingPairNode));
  end

  UP = 1 - DelProb;                     % probability of being used
  for nn = TertiaryFreeNode:n,
    UP = UP * 0.99;                     % reduce usage probability
    Node(nn).Delete = 1 - UP;           % increasing deletion probs
  end

  s = 0.99;

  m = n - TertiaryFreeNode;                    % number of extra basepairs
  m = max(m,10);                               % at least 10
  m = min(m,20);                               % no more than 20

  t = 0.999999;                                % provisional usage probability

  while (s-s^(m+1))/(1-s) + (s^m)*(t-t^(m+1))/(1-t) > m,
    t = t * 0.999;
  end

  UP = t;

  LastcWWText = Node(n).Comment(5:end);

  for nn = 1:m,
    UP = UP * 0.9;                    % reduce usage probability

    n = n+1;  
    Node(n).type        = 'Basepair';        % node type
    Node(n).nextnode    = n+1;               % index of next node in tree
    Node(n).LeftLetter  = 'C';
    Node(n).RightLetter = 'G';
    Node(n).Edge        = 1;
    Node(n).Delete      = 1-UP;        % deletion prob
    Node(n).lpar        = [.01*ones(16,1); 0]; % left insertion param
    Node(n).rpar        = [.01*ones(16,1); 0]; % right insertion param
    Node(n).LeftIndex   = a;
    Node(n).RightIndex  = B;

    Node(n).SubsProb=pIsoScore(1,'C','G',method,ExemplarIDI,ExemplarFreq);

    L = Node(n).lpar(1,1);
    R = Node(n).rpar(1,1);
    X = 0:10;     
    Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
                * sum(R.^X*exp(-R)./factorial(X));
    Node(n).leftLengthDist  = subPoisson(Node(n).lpar(1,1));
    Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
    Node(n).rightLengthDist = subPoisson(Node(n).rpar(1,1)); 
    Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

    Node(n).Comment = [' // Extra basepair node CG cWW'];

    if Verbose > 0,
      fprintf('%3d Extra basepair after %s\n', n, LastcWWText);
    end
  end

  for nn = TertiaryFreeNode:n,
    if Verbose > 1,
      fprintf('Node %3d has deletion probability %8.6f\n', nn, Node(nn).Delete);
    end
  end

  for nn = 1:length(ClosingPairNode),
    n = n+1;
    Node(n) = ClosingPairNode(nn);
    Node(n).nextnode = n + 1;
  end
