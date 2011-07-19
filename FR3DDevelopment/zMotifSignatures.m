% zMotifSignatures(Motif,Strands,Rotation,Type) determines a motif signature for the motif whose FR3D query information is encoded in Motif.  Strands is the number of strands in the motif, Rotation tells which strand to start with, Type is 0 for basepairs only, 1 for basepairs and stacks

function [Sig,AllSig] = zMotifSignatures(Motif,Strands,Rotation,Type)

if nargin < 4,
  Type = 0;                                       % basepairs only
end

if nargin < 3,
  Rotation = 1;
end

if nargin < 2,
  Strands = 2;
end

if ~isfield(Motif,'EdgeNums'),
  Motif.EdgeNums = Motif.Edge;
end

N = length(Motif.EdgeNums(1,:));

Sig = '';

% ----------------------------------------------- % determine strand boundaries

if Strands == 1,
  S = N-1;
elseif Strands == 2,


  if any(Motif.EdgeNums{1,N} == 1),          % flanking cWW pair
    i = N;
    CWW = 0;
    while CWW == 0 && i > 0,
      for j = (i+1):N,
        if any(Motif.EdgeNums{i,j} == 1),
          CWW = 1;
        end
      end
      if CWW == 0,
        i = i - 1;
      end
    end
                                             % probe for more in this strand
    while Motif.DifferenceSign(i,i+1) == -1 && i < N-1,
      i = i + 1;
    end
    S = i;
  else
    S = N;                                   % go through all rows
  end

  if Rotation == 2;
    p = [(S+1):N 1:S];                       % re-order strands
    Motif.EdgeNums = Motif.EdgeNums(p,p);
    Motif.Edges    = Motif.Edges(p,p);
    S = N-S;
  end

  S
  Motif.Edges

end

for i = 1:S,
  Sig = [Sig '-'];                                % separator between nts
  for j = N:-1:i,
    Sig = [Sig ','];                                % separator between nts
    if length(Motif.EdgeNums{i,j}) >= 1,           % one or more here
      e = Motif.EdgeNums{i,j};                  % edge codes
      if ~isempty(strfind(lower(Motif.Edges{i,j}),'pair')) ...
        || ~isempty(strfind(lower(Motif.Edges{j,i}),'pair')),
        Sig = [Sig 'pair'];
      else
        for k = 1:length(e),
          if Type == 1,
            Sig = [Sig zEdgeText(e(k))];           % add pairs and stacks
            if length(e) > k ,
              Sig = [Sig '/'];                  % multiple pairs here
            end
          elseif abs(e(k)) < 13,                   % only basepairs
            Sig = [Sig zEdgeText(e(k))];
            if length(e) > k ,
              Sig = [Sig '/'];                  % multiple pairs here
            end
          end
        end
      end
    end
  end
end

Sig = strrep(Sig,' ','');
for i = 1:N,
  Sig = strrep(Sig,',,',',');
end
Sig = strrep(Sig,',-','-');
Sig = strrep(Sig,'-,','-');
Sig = strrep(Sig,'cWW/cWW','cWW');
Sig = strrep(Sig,'tWW/tWW','tWW');
Sig = strrep(Sig,'cHH/cHH','cHH');
Sig = strrep(Sig,'tHH/tHH','tHH');
Sig = strrep(Sig,'cSS/cSS','cSS');
Sig = strrep(Sig,'tSS/tSS','tSS');
Sig = strrep(Sig,'s33/s33','s33');
Sig = strrep(Sig,'s55/s55','s55');

while Sig(1) == '-' && length(Sig) > 1,
  Sig = Sig(2:end);
end

while Sig(end) == '-' && length(Sig) > 1,
  Sig = Sig(1:(end-1));
end

while Sig(1) == ',' && length(Sig) > 1,
  Sig = Sig(2:end);
end

while Sig(end) == ',' && length(Sig) > 1,
  Sig = Sig(1:(end-1));
end

AllSig{Rotation} = Sig;

if Strands > 1 && Rotation < Strands,
  AllSig{Rotation+1} = zMotifSignatures(Motif,Strands,Rotation+1,Type);
end

if Rotation == 1,
  for i = 1:Strands,
    AllSig{i}
  end
  pause
end

