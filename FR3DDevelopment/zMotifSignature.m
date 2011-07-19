% zMotifSignatures(Motif,Strands,Rotation,Type) determines a motif signature for the motif whose FR3D query information is encoded in Motif.  Strands is the number of strands in the motif, Rotation tells which strand to start with, Type is 0 for basepairs only, 1 for basepairs and near pairs, 2 for pairs and stacks

function [Sig,AllSig] = zMotifSignature(Edge,Strands,Rotation,Type)

if nargin < 4,
  Type = 0;                                       % basepairs only
end

if nargin < 3,
  Rotation = 1;
end

if nargin < 2,
  Strands = 2;
end

N = length(Edge(1,:));                            % number of nucleotides

E = abs(fix(Edge));

if Type == 0,
  E = E .* (E < 14);                              % basepairs only
elseif Type == 1,
  E = E .* (E < 14) + E .* (E > 100) .* (E < 114); % basepairs and near pairs
elseif Type == 2,
  E = E .* (E < 24);                              % basepairs and stacks
end

Sig = '';                                         % signature starts blank

% ----------------------------------------------- % determine strand boundaries

if Strands == 1,
  S = N-1;
elseif Strands == 2,

  if E(1,N) == 1,                            % flanking cWW pair
    i = N;
    CWW = 0;
    while CWW == 0 && i > 0,
      for j = (i+1):N,
        if E(i,j) == 1,
          CWW = 1;
        end
      end
      if CWW == 0,
        i = i - 1;
      end
    end
                                             % probe for more in this strand
%    while Motif.DifferenceSign(i,i+1) == -1 && i < N-1,
%      i = i + 1;
%    end
    S = i;
  else
    S = N;                                   % go through all rows
  end

  if Rotation == 2;
    p = [(S+1):N 1:S];                       % re-order strands
    Edge = Edge(p,p);
    E    = E(p,p);
    S = N-S;
  end

end

for i = 1:N-1,                                % loop over rows of Edge
  Sig = [Sig '-'];                            % separator between nts
  for j = (i+1):N,                            % loop over columns
    Sig = [Sig ','];                          % separator between nts
    if E(i,j) > 0,                            % one or more here
      e = Edge(i,j);                          % edge codes
      for k = 1:length(e),
        Sig = [Sig zEdgeText(e(k))];          % add pairs and stacks
        if length(e) > k ,
          Sig = [Sig '/'];                    % multiple pairs here
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
Sig = strrep(Sig,'ncWW/ncWW','ncWW');
Sig = strrep(Sig,'ntWW/ntWW','ntWW');
Sig = strrep(Sig,'ncHH/ncHH','ncHH');
Sig = strrep(Sig,'ntHH/ntHH','ntHH');
Sig = strrep(Sig,'ncSS/ncSS','ncSS');
Sig = strrep(Sig,'ntSS/ntSS','ntSS');
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


Sig = [sprintf('%2d',N) '_' Sig];

Sig = strrep(Sig,' ','0');


AllSig{Rotation} = Sig;

if Strands > 1 && Rotation < Strands,
  AllSig{Rotation+1} = zMotifSignature(Edge,Strands,Rotation+1,Type);
end

% Anton 5/31/2011 Turned this off temporarily to make a cleaner log file
% if Rotation == 1,
%   for i = 1:Strands,
%     AllSig{i}
%   end
% end
% Anton 5/31/2011 Turned this off temporarily to make a cleaner log file

