% pConsensusPairSubstitution(a,b,f,File,F,L,Search) looks at the letter pairs corresponding to nucleotides a and b of the query motif in Search, uses interactions in F as the consensus, 

function [Score] = pConsensusPairSubstitution(a,b,f,File,F,L,Search,Param)

method = 2;             % default method for assigning pair subst probs

if length(Param) > 1,
  method  = Param(2);
end

N = length(Search.Candidates(1,:)) - 1;

Verbose = Param(1);

Verbose = 1;

load PairExemplars

Score = zeros(4,4);                      % ready to sum IsoScores for this pair 

cl = fix(F.Edge(a,b));                   % actual pair or stack between a and b

if any(cl == [-1 -2 -7 -8 -22 -23 -101 -102 -107 -108 -122 -123]),
  cl = -cl;                                   % symmetric pair, use positive
end

% -------------------------------------- collect interactions across instances

for c = 1:L,                                  % loop through candidates
  e = [];                                     % record observed interaction
  cp= [];
  for c = 1:L,                                % run through candidates
    i = Search.Candidates(c,a);               % index of first nucleotide
    j = Search.Candidates(c,b);               % index of second nucleotide
    e = [e Search.File(f(c)).Edge(i,j)];      % append observed interaction
    cp= [cp Search.File(f(c)).Coplanar(i,j)]; % append coplanar information
  end

  e = fix(e);                                 % round subcategories

  for d = 1:length(e),
    if any(e(d) == [-1 -2 -7 -8 -101 -102 -107 -108]),% don't distinguish sign
      e(d) = -e(d);
    end
  end
end

% -------------------------------------- if consensus on type of interaction

if sum(e == cl) >= 0.5*L && abs(cl)>0 && abs(cl) < 30 ,   
                                              % more than half in the categ
  e = ones(size(e)) * cl;                     % use this one every time
%  fprintf('pConsensusPairSubstitution: Consensus is %4s\n', zEdgeText(cl));
else
%  fprintf('pConsensusPairSubstitution: No consensus; near pairs will be used\n');  
end

fprintf('pConsensusPairSubstitution: Consensus is %4s\n', zEdgeText(cl));

% --------------------------------- Accumulate 4x4 matrices from instances

for c = 1:L,                            % loop through candidates
  i   = Search.Candidates(c,a);         % index 1 of pair in candidate
  j   = Search.Candidates(c,b);         % index 2 of pair in candidate
  NT1 = File(f(c)).NT(i);               % retrieve the first nucleotide
  NT2 = File(f(c)).NT(j);               % retrieve the second nucleotide

  if Verbose > 0,
    fprintf('pConsensusPairSubstitution: File %4s has %s%4s and %s%4s making %4s', File(f(c)).Filename, NT1.Base, NT1.Number, NT2.Base, NT2.Number, zEdgeText(File(f(c)).Edge(i,j)));
  end

  g = fix(File(f(c)).Edge(i,j));

  if any(g == [-1 -2 -7 -8 -22 -23 -101 -102 -107 -108 -122 -123]),
    g = -g;                                  % symmetric pair, use positive
  end

  if g == cl || g == sign(cl)*(100+abs(cl)),
    newScore = pIsoScore(e(c),NT1.Base,NT2.Base,method,ExemplarIDI,ExemplarFreq);
                                             % use consensus edge
    Score = Score + newScore;
    fprintf(', using it\n');
  else
    fprintf(', not using\n');
  end

%newScore(1:7)
end

Score = Score / L;                          % average *** need better idea!

%Score
