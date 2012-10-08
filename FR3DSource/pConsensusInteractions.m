% pConsensusInteractions(Search) determines the consensus interactions in the candidates in Search

function [Edge,BPh,BR,Search] = pConsensusInteractions(Search)

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

f = Search.Candidates(:,N+1);           % file numbers of motifs

BPh = zeros(N,N);
BR  = zeros(N,N);

% ----- Calculate coplanar measure for each File in Search

if ~isfield(Search.File(1),'Coplanar'),

  disp('Adding coplanar values')

  clear NewFile
  for ff = 1:length(Search.File),
    F = Search.File(ff);
    if ~isempty(F.NT),
      F.Coplanar = sparse(F.NumNT,F.NumNT);
      NewFile(ff) = F;
    end
  end
  Search.File = NewFile;

  for c = 1:length(Search.Candidates(:,1)),
    ff = Search.Candidates(c,N+1);
    i = Search.Candidates(c,1:N);
    for a = 1:length(i),
      for b = (a+1):length(i),
        if Search.File(ff).Edge(i(a),i(b)) ~= 0,
          NT1 = Search.File(ff).NT(i(a));
          NT2 = Search.File(ff).NT(i(b));
          Pair = zAnalyzePair(NT1,NT2);
          Search.File(ff).Coplanar(i(a),i(b)) = Pair.Coplanar;
          Search.File(ff).Coplanar(i(b),i(a)) = Pair.Coplanar;
        end
      end 
    end
  end
end

% -------------------------------------------- find consensus pairs and stacks

Edge = sparse(N,N);                             % place to store consensus

for a = 1:N,                                    % first NT of possible pair
  for b = (a+1):N,                              % second NT of possible pair
    e = [];                                     % record observed edges
    cp= [];                                     % coplanar, for basepairs
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).Edge(i,j)];      % append observed interaction
      
      % modified by Anton to deal with old motif atlas releases 12/6/2011
      if isfield(Search.File(f(c)),'LooseCoplanar')
          cp= [cp Search.File(f(c)).LooseCoplanar(i,j)]; % append coplanar information
      else
          cp= [cp 0];
      end
     % modified by Anton to deal with old motif atlas releases 12/6/2011

    end

    e = fix(e);                                 % round subcategories

    for d = 1:length(e),                        % loop through instances
      if any(e(d) == [-1 -2 -7 -8 -22 -23 -101 -102 -107 -108 -122 -123]),
        e(d) = -e(d);                           % don't distinguish sign here
      end
    end

    ae = (abs(e)>0) .* (abs(e) < 30) .* e;      % FR3D basepairs and stacks

    anc= (abs(e)>100) .* (abs(e) < 120) .* (cp > 0) .* e; % npair but coplanar

    % ------------------------------------------- determine the "best" inter

    if length(ae(ae~=0)) > 0,                   % if there is a true interact
      emode = mode(ae(ae~=0));                  % store the code for it
    elseif length(anc(anc~=0)) > 0,             % if there is a near coplanar
      emode = mode(anc(anc~=0));                % store the code for it
      emode = sign(emode)*(abs(emode)-100);     % corresponding true code
    else
      emode = NaN;
    end

%[a b]
%cp

%emode

    nemode = sign(emode) * (abs(emode) + 100);  % near category code
%nemode

    numeqemode = sum(emode == e);               % number equal to the mode
%numeqemode

    numeqnemode = sum(nemode == e);             % number equal to near mode

    numeqnemodec= sum(nemode == anc);           % coplanar and equal near mode
%numeqnemode

    if numeqemode >= max(min(L,2),L/8),         % enough true pairs
      Edge(a,b) = emode;                        % use the most common one
    elseif (numeqemode+numeqnemodec) > L/2,     % true and coplanar near
      Edge(a,b) = emode;                        % use the most common one
    elseif numeqemode >= L/3 && numeqnemode >= L/3, % some true, some near
      Edge(a,b) = emode;                        % use the most common one
    end

    Edge(b,a) = -Edge(a,b);                     % record twice

    % --------------------------------- old method, for comparison

    if any((abs(e)>0).*(abs(e)<30)),            % there was some bp or stack
      ee = e(find((abs(e) > 0) .* (abs(e) < 30))); % FR3D basepairs and stacks
      if ~isempty(ee),
        if mode(ee) ~= Edge(a,b),
%          fprintf('************ Consensus for (%d,%d) was %d now is %d\n', a,b, full(mode(ee)), full(Edge(a,b)));
        end
      end
    else                                        % only near interactions
      ee = e(find((abs(e) > 100) .* (cp > 0))); % near, coplanar pairs
      if ~isempty(ee),
        if mode(ee) ~= Edge(a,b),
%          fprintf('************ Consensus for (%d,%d) was %d now is %d\n', a,b, full(mode(ee)), full(Edge(a,b)));
        end
      end
    end
  end
end

% ----------------------------------------------- BPh and BR edges used

% BBEdge{1} = 

% ----------------------------------------------- consensus BPh interactions

for a = 1:N,                                    % first NT of possible BPh
  for b = 1:N,                                  % second NT of possible BPh
    e = [];                                     % record observed edges
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).BasePhosphate(i,j)];% append observed interaction
    end

    ae = (abs(e)>0) .* (abs(e) < 20) .* e;      % valid BPh interactions

    if length(find(ae)) > length(e)/2,          % more than half make BPh
      BPh(a,b) = 1;                             % indicator of conservation
%full(e)
%full(ae)
%BPh(a,b)
    end
  end
end

% ----------------------------------------------- consensus BR interactions

if ~isfield(Search.File(1),'BaseRibose'),
  disp('Adding base-ribose interactions')
  Search.File = zBaseRiboseInteractions(Search.File);
end

for a = 1:N,                                    % first NT of possible BR
  for b = 1:N,                                  % second NT of possible BR
    e = [];                                     % record observed edges
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).BaseRibose(i,j)];  % append observed interaction
    end

    ae = (abs(e)>0) .* (abs(e) < 20) .* e;      % valid BR interactions

    if length(find(ae)) > length(e)/2,          % more than half make BR
      BR(a,b) = mode(ae(ae > 0));

  % a more sophisticated method for determining the consensus is needed
  % count how many times each one occurs, giving maybe 0.5 for near
    end
  end
end

end
