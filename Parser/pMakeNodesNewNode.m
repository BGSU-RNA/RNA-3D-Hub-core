
Node(n).type      = '';             % node type
Node(n).nextnode  = 0;                   % index of next node in tree
Node(n).LeftIndex = 1;
Node(n).RightIndex= 1;
Node(n).Insertion = [];                    % make sure this field exists
Node(n).id        = '';                    % make sure this field exists
Node(n).leftLengthDist = [];
Node(n).leftLetterDist = [];
Node(n).rightLengthDist = [];
Node(n).rightLetterDist = [];
Node(n).LeftLetter = '';
Node(n).RightLetter = '';
Node(n).Comment = '';
Node(n).NumLoops = [];
if n == 1,
  Node(n).Edge = sparse(zeros(length(File.NT)));
  Node(n).Extensibility = zeros(1,length(File.NT));
end
Node(n).Delete = [];
Node(n).SubsProb = [];
Node(n).Z = [];
Node(n).MiddleIndex = [];
Node(n).P = [];
Node(n).PIns = [];
Node(n).subtype = [];
Node(n).Left = [];
Node(n).Right = [];
Node(n).InsertionComment = '';
Node(n).IBases = [];
Node(n).InteractionComment = '';
Node(n).JunctionDeletion = [];
Node(n).InterIndices = [];
