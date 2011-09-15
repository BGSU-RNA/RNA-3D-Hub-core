% set up basepair

n = n+1;  
Node(n).type        = 'Basepair';        % node type
Node(n).nextnode    = n+1;               % index of next node in tree
Node(n).LeftLetter  = File.NT(a).Base;
Node(n).RightLetter = File.NT(B).Base;
Node(n).Edge        = Interact{a}.Categ(1);
Node(n).Delete      = DelProb + (B-a)/1000000;        % deletion prob
Node(n).lpar        = [.01*ones(16,1); 0]; % left insertion param
Node(n).rpar        = [.01*ones(16,1); 0]; % right insertion param
Node(n).LeftIndex   = a;
Node(n).RightIndex  = B;
Node(n).leftLengthDist  = subPoisson(0.01);
Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
Node(n).rightLengthDist = subPoisson(0.01); 
Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

P = pIsoScore(File.Edge(a,B),Node(n).LeftLetter, ...
                         Node(n).RightLetter,method,ExemplarIDI,ExemplarFreq);

Node(1).Edge(a,B) = File.Edge(a,B);      % store this interaction

L = Node(n).lpar(1,1);
R = Node(n).rpar(1,1);

X = 0:10;                                % 0 to 10 insertions after a basepair

Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
          * sum(R.^X*exp(-R)./factorial(X));

Node(n).Comment = [' // Basepair ' File.NT(a).Base File.NT(a).Number ' - ' File.NT(B).Base File.NT(B).Number ' ' zEdgeText(File.Edge(a,B))];

if Verbose > 0,
  fprintf('%3d Basepair  %4s %4s %s%s %s\n',n, File.NT(a).Number, File.NT(B).Number,File.NT(a).Base,File.NT(B).Base,zEdgeText(File.Edge(a,B)));
end

% ----------------------------------------- Adjust subs probs for LR, BPh inter

if AdjustSubsForLR == 1 && (TertiaryFreeNode == 0 || Extension < 2),
  Node(n).SubsProb = pAdjustSubsProb(File,a,B,P,method);
else
  Node(n).SubsProb = P;
end

a = a + 1;                                % current base on the left
B = B - 1;                                % current base on the right

