%<div class="moz-text-flowed" style="font-family: -moz-fixed">


n=1;                                            % current node
  Node(n).type          = 'Initial';                % node type
  Node(n).nextnode      = n+1;                      % index of next node in tree
  Node(n).lpar          = 1.5;                        % left insertion parameter
  Node(n).rpar          = 1.5; 
  Node(n).Bl            = [1,3];
  Node(n).Br            = [120,122];
  Node(n).basehalfwidth = 6;

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;
  Node(n).Delete      = 0.1;
  Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [4];
  Node(n).Br        = [119];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;
  Node(n).Delete      = 0.1;
  Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [4];
  Node(n).Br        = [119];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.1;
  Node(n).lpar        = [.001*ones(16,1); 0];        % left insertion parameters  
  Node(n).rpar        = [.001*ones(16,1); 0];     % right insertion parameters  
  Node(n).Bl        = [4];
  Node(n).Br        = [119];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.012;
  Node(n).lpar        = [0.041*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.041*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [4];
  Node(n).Br        = [119];

  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.04*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.04*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [5];
  Node(n).Br        = [118];

  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.01;
  Node(n).lpar        = [.01*ones(16,1); 0];        % left insertion parameters  
  Node(n).rpar        = [.01*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [6];
  Node(n).Br        = [117];

  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.01;
  Node(n).lpar        = [.01*ones(16,1); 0];        % left insertion parameters  
  Node(n).rpar        = [.01*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [7];
  Node(n).Br        = [116];

  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.01;
  Node(n).lpar        = [.01*ones(16,1); 0];        % left insertion parameters  
  Node(n).rpar        = [.01*ones(16,1); 0];     % right insertion parameters  
  Node(n).Bl        = [8];
  Node(n).Br        = [115];

  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.01;
  Node(n).lpar        = [.005*ones(16,1); 0];        % left insertion parameters  
  Node(n).rpar        = [.005*ones(16,1); 0];     % right insertion parameters  
  Node(n).Bl        = [9];
  Node(n).Br        = [114];

  
  n=n+1;                                            % current node
  % a junction is here..it is anotated below the rest of the nodes


% 11-68 ---------------------------------------------------------

n=n+1;          % current node
  a=n;
  Node(n).type      = 'Initial';                % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).lpar      = 0.01;                     % left insertion parameter
  Node(n).rpar      = 0.01;                     % right insertion parameter  
  Node(n).Bl        = [16];
  Node(n).Br        = [64];
  
n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015;
  Node(n).lpar        = [0.016*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.016*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [16];
  Node(n).Br        = [64];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015;
  Node(n).lpar        = [0.015*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.015*ones(16,1); 0];     % right insertion parameters 
  Node(n).Bl        = [17];
  Node(n).Br        = [63];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'U';
  Node(n).RightLetter = 'A';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015;
  Node(n).lpar        = [0.014*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.014*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [18];
  Node(n).Br        = [62];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015;
  Node(n).lpar        = [0.013*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.013*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [19];
  Node(n).Br        = [61];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015; 
  Node(n).lpar        = [0.012*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.012*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl        = [20];
  Node(n).Br        = [60];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015; 
  Node(n).lpar        = [0.011*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.011*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [21];
  Node(n).Br          = [59];

n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1 4 5];                    % nucleotides to use on left
  Node(n).Left(2,:) = [1 3 4];     
  Node(n).Left(3,:) = [1 2 3];     
  Node(n).Left(4,:) = [1 5 6];     
  Node(n).LIP       = [0.033 0.033 0.033 0.9];  % probs for insertion possibs
  Node(n).Right(1,:)= [1];                      % nucleotides to use on right
  Node(n).RIP       = [1];                        % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 2];
  Node(n).Score(:,:,1) = pIsoScore(6.0,3,2);
  Node(n).IBases(2,:)  = [3 4];
  Node(n).Score(:,:,2) = pIsoScore(1.0,2,3);
  Node(n).Bl          = [22,27];
  Node(n).Br          = [58];

n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1 2];                    % nucleotides to use on left
  Node(n).Left(2,:) = [1 3];     
  Node(n).LIP       = [0.9 0.1];                % probs for insertion possibs
  Node(n).Right(1,:)= [5 4 2 1];                % nucleotides to use on right
  Node(n).Right(2,:)= [4 3 2 1];                % nucleotides to use on right
  Node(n).RIP       = [0.9 0.1];                % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 6];
  Node(n).Score(:,:,1) = pIsoScore(6.1,4,1);
  Node(n).IBases(2,:)  = [2 5];
  Node(n).Score(:,:,2) = pIsoScore(5.0,2,1);
  Node(n).IBases(3,:)  = [1 4];
  Node(n).Score(:,:,3) = pIsoScore(1,4,1);
  Node(n).IBases(4,:)  = [2 3];
  Node(n).Score(:,:,4) = pIsoScore(1,2,3);
  Node(n).Bl          = [28,29];
  Node(n).Br          = [53,57];

n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1];                      % nucleotides to use on left
  Node(n).LIP       = [1];                      % probs for insertion possibs
  Node(n).Right(1,:)= [3 2 1];                  % nucleotides to use on right
  Node(n).Right(2,:)= [4 3 1];                  % nucleotides to use on right
  Node(n).RIP       = [0.9 0.1];                % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 4];
  Node(n).Score(:,:,1) = pIsoScore(5.0,2,1);
  Node(n).IBases(2,:)  = [3 4];
  Node(n).Score(:,:,2) = pIsoScore(9.0,1,1);
  Node(n).IBases(3,:)  = [1 2];
  Node(n).Score(:,:,3) = pIsoScore(1.0,2,3);
  Node(n).Bl          = [30];
  Node(n).Br          = [50,52];

n=n+1;                                            % current node
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.005;
  Node(n).lpar        = [0.01*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.01*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [31];
  Node(n).Br          = [49];

n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1 2 3 4];                    % nucleotides to use on left
  Node(n).Left(2,:) = [1 2 4 5];     
  Node(n).LIP       = [0.1 0.9];                % probs for insertion possibs
  Node(n).Right(1,:)= [2 1];                      % nucleotides to use on right
  Node(n).RIP       = 1;                        % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 6];
  Node(n).Score(:,:,1) = pIsoScore(1.0,3,2);
  Node(n).IBases(1,:)  = [2 5];
  Node(n).Score(:,:,1) = pIsoScore(1.0,4,1);
  Node(n).IBases(2,:)  = [4 5];
  Node(n).Score(:,:,2) = pIsoScore(-5.0,2,1);
  Node(n).IBases(2,:)  = [3 6];
  Node(n).Score(:,:,2) = pIsoScore(-5.0,2,2);
  Node(n).Bl          = [32,36];
  Node(n).Br          = [47,48];

n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1 2 3];                  % nucleotides to use on left
  Node(n).Left(2,:) = [1 3 4];                  % nucleotides to use on left
  Node(n).Left(3,:) = [1 2 4];                  % nucleotides to use on left
  Node(n).LIP       = [0.8 0.1 0.1];            % probs for insertion possibs
  Node(n).Right(1,:)= [4 3 2 1];                % nucleotides to use on right
  Node(n).RIP       = 1;                        % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 4];
  Node(n).Score(:,:,1) = pIsoScore(1.0,2,3);
  Node(n).IBases(2,:)  = [2 6];
  Node(n).Score(:,:,2) = pIsoScore(5.0,1,1);
  Node(n).IBases(3,:)  = [3 5];
  Node(n).Score(:,:,3) = pIsoScore(10.0,4,1);
  Node(n).IBases(4,:)  = [4 7];
  Node(n).Score(:,:,4) = pIsoScore(6.0,3,2);
  Node(n).Bl          = [37,39];
  Node(n).Br          = [43,46];
  
n=n+1;                                            % current node
  Node(n).type      = 'Hairpin';                % node type
  Node(n).subtype   = 'None';                   % nothing special
  Node(n).nextnode  = Inf;                      % index of next node in tree
  Node(n).P         = ones(17,1);               % state to state transitions
  Node(n).PIns      = 1;                        % when no previous state
  Node(n).Bm          = [40,42];


%   % 69-112 ---------------------------------------------------------------------

 n=n+1;    % current node
 b=n;
  Node(n).type        = 'Initial';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).lpar        = 0.005;                     % left insertion parameter
  Node(n).rpar        = 0.005;                     % right insertion parameter
  Node(n).Bl          = [69];
  Node(n).Br          = [112];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'U';
  Node(n).RightLetter = 'U';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [69];
  Node(n).Br          = [112];

 n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'U';
  Node(n).RightLetter = 'U';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters  
  Node(n).Bl          = [70];
  Node(n).Br          = [111];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.012;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [71];
  Node(n).Br          = [110];

 n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [72];
  Node(n).Br          = [109];

 n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [73];
  Node(n).Br          = [108];

 n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [74];
  Node(n).Br          = [107];

  n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.1*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.1*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [75];
  Node(n).Br          = [106];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'A';
  Node(n).Inter       = 10;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [76];
  Node(n).Br          = [105];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'A';
  Node(n).RightLetter = 'A';
  Node(n).Inter       = 8;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [77];
  Node(n).Br          = [104];
  
 n=n+1;                                            % current node
  Node(n).type      = 'Motif';                  % node type
  Node(n).nextnode  = n+1;                      % index of next node in tree
  Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
  Node(n).PIns	    = [0.05 0.95];              % when no previous state
  Node(n).Left(1,:) = [1 2];                    % nucleotides to use on left
  Node(n).LIP       = [1];                      % probs for insertion possibs
  Node(n).Right(1,:)= [1];                      % nucleotides to use on right
  Node(n).RIP       = [1];                      % probs for insertion possibs
  Node(n).IBases(1,:)  = [1 2];
  Node(n).Score(:,:,1) = pIsoScore(-9,3,4);
  Node(n).IBases(2,:)  = [2 3];
  Node(n).Score(:,:,2) = pIsoScore(-4,4,1);
  Node(n).Bl          = [78,79];
  Node(n).Br          = [103];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'A';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 10;  
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [80];
  Node(n).Br          = [102];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'C';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;
  Node(n).Delete      = 0.011;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [81];
  Node(n).Br          = [101];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'U';
  Node(n).RightLetter = 'G';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.012;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [82];
  Node(n).Br          = [100];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'U';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.013;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [83];
  Node(n).Br          = [99];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.014;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [84];
  Node(n).Br          = [98];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'A';
  Node(n).RightLetter = 'U';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.015;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [85];
  Node(n).Br          = [97];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;  
  Node(n).Delete      = 0.016;
  Node(n).lpar        = [1*ones(16,1); 0];        % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [86,87];
  Node(n).Br          = [96];

n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1; 
  Node(n).Delete      = 0.017;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [88];
  Node(n).Br          = [95];
  
n=n+1;
  Node(n).type        = 'Basepair';                % node type
  Node(n).nextnode    = n+1;                      % index of next node in tree
  Node(n).LeftLetter  = 'G';
  Node(n).RightLetter = 'C';
  Node(n).Inter       = 1;
  Node(n).Delete      = 0.018;
  Node(n).lpar        = [0.005*ones(16,1); 0];     % left insertion parameters
  Node(n).rpar        = [0.005*ones(16,1); 0];     % right insertion parameters
  Node(n).Bl          = [89];
  Node(n).Br          = [94];

n=n+1;                                           % current node
  Node(n).type      = 'Hairpin';                % node type
  Node(n).subtype   = 'GNRA';                   % GNRA loop
  Node(n).nextnode  = Inf;                      % index of next node in tree
  Node(n).P         = ones(17,1);               % state to state transitions
  Node(n).PIns      = 1;                        % when no previous state
  Node(n).Bm          = [90,93];


  N=n;
  
  n=a-1;
   Node(n).type      = 'JunctionMotif';               % node type
   Node(n).nextnode  =  [a b];                   % index of next node in tree
   Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                % state to state transitions
   Node(n).PIns	        = [0.05 0.95];                                % when no previous state
   Node(n).Left(1,:)    = [1 2 4 5 6];                                   % nucleotides to use on left
   Node(n).Left(2,:)    = [1 2 3 4 5];                                   % nucleotides to use on left
   Node(n).LIP          = [.98 .02];    % probs for insertion possibs
   Node(n).Middle(1,:)  = [1 2 3 4];                                 % nucleotides to use in middle
   Node(n).MIP          = [1];                                          % probs for insertion possibs
   Node(n).Right(1,:)   = [1];                                       % nucleotides to use on right
   Node(n).RIP          = [1];                                         % probs for insertion possibs
   Node(n).Letters      ='CAAGCAGCGC';
   Node(n).IBases(1,:)  = [1 3];
   Node(n).Inter(1)     = 6;
   Node(n).IBases(2,:)  = [2 9];
   Node(n).Inter(2)     = 12;
   Node(n).IBases(3,:)  = [4 8];
   Node(n).Inter(3)     = 1;
   Node(n).IBases(4,:)  = [5 7];
   Node(n).Inter(4)     = 1;
   Node(n).IBases(5,:)  = [6 8];
   Node(n).Inter(5)     = -6;
   Node(n).IBases(6,:)  = [7 10];
   Node(n).Inter(6)     = 4.1;
   Node(n).Bl        = [10,15];
   Node(n).Bm        = [65,68];
   Node(n).Br        = [113];
%</div>