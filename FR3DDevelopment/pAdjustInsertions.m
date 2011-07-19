% Adjust insertion proabilities and letter distributions
% For use with basepair nodes and initial nodes following cluster nodes

% ----------------------------- tally insertions on the left
    inscount = [];
    letter = Prior;                      % record which bases occur
    if n < length(Node)-1,                    % no insertions after last pair
      aa = min(Node(n+1).LeftIndex);          % next interacting base in the model
      for c = 1:L,
        inscount(c) = abs(Search.Candidates(c,aa) - Search.Candidates(c,a)) - 1;
      end

%disp('pMakeMotifModelFromSSF: Left insertion counts:')
%inscount

      lld = zeros(1,max(inscount)+2);           % Dirichlet distribution
      for c = 1:L,
        lld(inscount(c)+1) = lld(inscount(c)+1) + 1;
        ff = Search.Candidates(c,N+1);          % file number
        d = Search.Candidates(c,a);             % index of interacting base
        if inscount(c) > 0,
           for i = 1:inscount(c),
             insb = Search.File(ff).NT(d+i).Code;  % A=1, C=2, G=3, U=4
             if ~isempty(insb)
               letter(insb) = letter(insb) + 1;
             end
           end
        end
      end

%disp('pMakeMotifModelFromSSF: Left insertion tallies:')
%lld
      numz = sum(lld==0);
      lld(lld==0) = (.01*sum(lld))/(1-.01*numz);
      Node(n).leftLengthDist = lld / sum(lld);    % normalize
      Node(n).leftLetterDist = letter / sum(letter);  % normalize 

%     bb = max(Node(n+1).RightIndex);
%     [rld,letter] = pInsertionConcenses(Search,Node,n,bb,b,Prior);
%     Node(n).rightLengthDist = lld;
%     Node(n).rightLetterDist = letter;

      % ----------------------------- tally insertions on the right
      inscount = [];
      letter = Prior;
      bb = max(Node(n+1).RightIndex);      % next interacting base in the model
      if bb==b
        inscount(1:L) = 0;
      else
        for c = 1:L,
          inscount(c) = abs(Search.Candidates(c,b) - Search.Candidates(c,bb)) - 1;
        end
      end
    
%disp('pMakeMotifModelFromSSF: Right insertion counts:')
%inscount
      rld = zeros(1,max(inscount)+2);  % Dirichlet distribution
      for c = 1:L,
        rld(inscount(c)+1) = rld(inscount(c)+1) + 1;
        ff = Search.Candidates(c,N+1);          % file number
        d = Search.Candidates(c,b);            % index of interacting base
        if inscount(c) > 0,
           for i = 1:inscount(c),
             insb = Search.File(ff).NT(d-i).Code; % A=1, C=2, G=3, U=4
             if ~isempty(insb)
               letter(insb) = letter(insb) + 1;
             end
           end
        end
      end
%disp('pMakeMotifModelFromSSF: Right insertion tallies:')
%rld
      numz = sum(rld==0);
      rld(rld==0) = (.01*sum(rld))/(1-.01*numz);
      rld(rld==0) = sum(rld)*.01;
      Node(n).rightLengthDist = rld / sum(rld);    % normalize
      Node(n).rightLetterDist = letter / sum(letter);  % normalize 
    end