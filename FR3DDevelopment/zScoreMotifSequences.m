% zScoreMotifSequences scores the sequences in codes using the consensus interactions found in Edge and BPh, using the given basepair method.  Full Search information is necessary to match the geometry of the stacks in the candidates.

function [scores,Name,subscores,SSNames] = zScoreMotifSequences(codes,Search,Edge,BPh,method,Stacks,BR)

Verbose = 1;

Name{1} = ['Basepairs using method ' num2str(method)];
Name{2} = ['Basepairs using method ' num2str(method) ' and base-phosphates'];
Name{3} = ['Basepairs using ' num2str(method) ', base-phosphates, and stacking'];
Name{4} = ['Basepairs using ' num2str(method) ', base-phosphates, base-ribose, and stacking'];
Name{5} = ['Exp of negative Hamming distance minus uniform random'];

[N,M] = size(Edge);
L     = length(Search.Candidates(:,1));         % number of candidates
[a,b] = size(codes);

for k = 1:L,
  f = Search.Candidates(k,N+1);                 % file number
  for n = 1:N,
    CandCodes(k,n) = Search.File(f).NT(Search.Candidates(k,n)).Code;
  end
end

% -------------------------------- score all possible sequences for basepairs

m = 1;
scores = ones(a,1);                    % place to store all scores
subscores = ones(a,1);
SSCount = 1;                           % which subscore we're on

if 10 > 1,

for p = 1:N,                                    % number of nucleotides
  for q = (p+1):N,                              % loop through interactions
    if abs(Edge(p,q)) > 0 && abs(Edge(p,q)) < 13,       % basepair
      IS = zeros(4,4);
      for k = 1:L,                              % loop through instances
        IS = IS + pIsoScore(Edge(p,q),CandCodes(k,p),CandCodes(k,q),method);
      end
      IS = IS / sum(sum(IS));

      if Verbose > 0,
        fprintf('Nucleotides %d and %d make %s, giving this substitution matrix:\n', p, q, zEdgeText(Edge(p,q)));
        IS
      end

      for s = 1:a,                              % loop through all sequences
        scores(s,m) = scores(s,m) * IS(codes(s,p),codes(s,q));
        subscores(s,SSCount) = IS(codes(s,p),codes(s,q));
      end

      SSNames{SSCount} = ['Basepair NT ' num2str(p) ' and ' num2str(q)];

      SSCount = SSCount + 1;

    end
  end
end

end

% -------------------------------- score all possible sequences for BPh

m = 2;
scores(:,m) = scores(:,m-1);                  % place to store all scores

if 10 > 1,

for p = 1:N,                                      % loop through bases
  for q = 1:N,
    if p~=q && abs(BPh(p,q)) > 0 && abs(BPh(p,q)) < 20,  % BPh consensus here
      P = zeros(1,4);                             % BPh substitution probs
      for k = 1:L,                                % loop through instances
        f = Search.Candidates(k,N+1);
        B = Search.File(f).BasePhosphate;
        x = Search.Candidates(k,p);
        y = Search.Candidates(k,q);
        [D,Q] = zBasePhosphateGeometryOld(mod(B(x,y),100));

        if ~isempty(Q),
          P = P + Q;
        end
      end

      P = P / sum(P);                             % normalize

      if Verbose > 0,
          fprintf('%d %d We get this BPh data\n');

[p q]
P
      end

      for s = 1:a,                      % loop through all sequences
        scores(s,m) = scores(s,m) * P(codes(s,p));
        subscores(s,SSCount) = P(codes(s,p));
      end

      SSNames{SSCount} = ['BPh NT ' num2str(p) ' and ' num2str(q)];

      SSCount = SSCount + 1;

    end
  end
end

end

% -------------------------------- score all possible sequences for stacking

m = 3;
scores(:,m) = scores(:,m-1);                  % place to store all scores

for p = 1:N,                                      % loop through bases
  for q = (p+1):N,                                
    if abs(Edge(p,q)) > 20 && abs(Edge(p,q)) < 24,  % base stack

      switch fix(Edge(p,q)),
      case 21,
        S = Stacks{1};                          % s35
      case {22, -22},
        S = Stacks{2};                          % s33
      case {23, -23},
        S = Stacks{3};                          % s55
      case -21,
        S = Stacks{1};                          % s53
        S.Candidates = S.Candidates(:,[2 1 3]);   % swap nucleotide order
      otherwise,
        disp('What case is this?');
        fix(Edge(p,q))
      end

% S.SaveName

      NS = length(S.Candidates(:,1));           % number of stacks
      for c = 1:NS,
        f = S.Candidates(c,3);                  % file number
        SC(c,1) = S.File(f).NT(S.Candidates(c,1)).Code;  % code of first base
        SC(c,2) = S.File(f).NT(S.Candidates(c,2)).Code;  % code of second base
      end

      H = zeros(4,4);

      for k = 1:L,                              % loop through instances
        f = Search.Candidates(k,N+1);           % file number

        clear Model

        Model.NT(1) = Search.File(f).NT(Search.Candidates(k,p));
        Model.NT(2) = Search.File(f).NT(Search.Candidates(k,q));

% zDisplayNT(Model)

        [Discrepancy, Candidates, i] = xRankCandidatesIDI(S.File,Model,S.Candidates,0);

        NewSC = SC(i,:);                        % re-order candidate codes

        MinIDI = 40*ones(4,4);
        
        for c = NS:-1:1,
          MinIDI(NewSC(c,1),NewSC(c,2)) = Discrepancy(c);
        end

        G = 1 ./ (1+MinIDI.^2);
        G = G / sum(sum(G));


        if Verbose > 0,
          fprintf('%d %d We get this min IDI matrix for %s %s%s %s%s %s %s\n',p,q,Search.File(f).Filename,Model.NT(1).Base,Model.NT(1).Number,Model.NT(2).Base,Model.NT(2).Number,zEdgeText(fix(Edge(p,q))));
 
MinIDI

          fprintf('We turn it into this matrix of scores:\n');
          G
        end

        H = H + G;

      end

      H = H / sum(sum(H))

      if Verbose > 0,
        fprintf('Overall matrix for the stack between %d and %d is\n',p,q);
        H
      end

      for s = 1:a,                      % loop through all sequences
        scores(s,m) = scores(s,m) * H(codes(s,p),codes(s,q));
        subscores(s,SSCount) = H(codes(s,p),codes(s,q));
      end

      SSNames{SSCount} = ['Stack NT ' num2str(p) ' and ' num2str(q)];

      SSCount = SSCount + 1;

    end
  end
end

% -------------------------------- score all possible sequences for Base-ribose

m = 4;
scores(:,m) = scores(:,m-1);                  % place to store all scores

if 10 > 1,

for p = 1:N,                                      % loop through bases
  for q = 1:N,
    if p~=q && abs(BR(p,q)) > 0 && abs(BR(p,q)) < 20,  % BR consensus here
      P = zeros(1,4);                             % BR substitution probs
      for k = 1:L,                                % loop through instances
        f = Search.Candidates(k,N+1);
        B = Search.File(f).BaseRibose;
        x = Search.Candidates(k,p);
        y = Search.Candidates(k,q);
        [D,Q] = zBasePhosphateGeometryOld(mod(B(x,y),100));

        if ~isempty(Q),
          P = P + Q;
        end
      end

      P = P / sum(P);                             % normalize

      if Verbose > 0,
          fprintf('%d %d We get this base-ribose data\n');

[p q]
P
      end

      for s = 1:a,                      % loop through all sequences
        scores(s,3) = scores(s,3) * P(codes(s,p));
        subscores(s,SSCount) = P(codes(s,p));
      end

      SSNames{SSCount} = ['BR NT ' num2str(p) ' and ' num2str(q)];

      SSCount = SSCount + 1;

    end
  end
end

end

% ------------------------------------------- normalize scores in each column

for k = 1:length(scores(1,:)),
  scores(:,k) = scores(:,k) / sum(scores(:,k));
end

end
