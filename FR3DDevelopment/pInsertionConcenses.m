% Changes the insertion probabilities of Node(n) between columns a and aa
% based on the sequences in Search 

% a <= aa

function [ld,letter] = pInsertionConcenses(Search,Node,n,a,aa,Prior)
    [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
    N = N - 1;                              % number of nucleotides

    % ----------------------------- tally insertions on the left
    inscount = [];
    letter = Prior;                      % record which bases occur
    if n < length(Node),
        inscount(1:L) = 0;
      if aa==a,
            inscount(1:L) = 0;
          for c = 1:L,
            inscount(c) = Search.Candidates(c,aa) - Search.Candidates(c,a) - 1;
          end
      end
    end

disp('pInsertionConcenses: Insertion counts:')
inscount

    ld = zeros(1,max(inscount)+2);           % Dirichlet distribution
    for c = 1:L,
      ld(inscount(c)+1) = ld(inscount(c)+1) + 1;
      ff = Search.Candidates(c,N+1);          % file number
      d = Search.Candidates(c,a);            % index of interacting base
      if inscount(c) > 0,
         for i = 1:inscount(c),
           insb = Search.File(ff).NT(d+i).Code;  % A=1, C=2, G=3, U=4
           if ~isempty(insb)
             letter(insb) = letter(insb) + 1;
           end
         end
      end
    end

disp('pInsertionConcenses: Insertion tallies:')
ld
    numz = sum(ld==0);
    ld(ld==0) = (.01*sum(ld))/(1-.01*numz);
    ld = ld / sum(ld);    % normalize
    letter = letter / sum(letter);  % normalize