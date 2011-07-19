% pSequenceConstraints puts limits on what each node could generate

function Sequence = pSequenceConstraints(Node,Sequence)

M = length(Sequence);                    % number of sequences
N = length(Node);                        % number of nodes

for k=1:M;                               % loop through all M sequences
    
  L = length(Sequence(k).X);             % number of letters in sequence k

  for n = 1:N                            % loop through nodes
    if ~isempty(Node(n).LeftIndex) % we know what Node n generates in 3D struct
      Sequence(k).Constraint(n).LeftIndex = ...
      Sequence(k).ColumntoIndex(Sequence(1).IndextoColumn(Node(n).LeftIndex));
    else
      Node(n).LeftIndex = 1;
    end

    if ~isempty(Node(n).MiddleIndex)
      Sequence(k).Constraint(n).MiddleIndex = ...
      Sequence(k).ColumntoIndex(Sequence(1).IndextoColumn(Node(n).MiddleIndex));
    else
      Node(n).MiddleIndex = 1;
    end

    if ~isempty(Node(n).RightIndex)
      Sequence(k).Constraint(n).RightIndex = ...
      Sequence(k).ColumntoIndex(Sequence(1).IndextoColumn(Node(n).RightIndex));
    else
      Node(n).RightIndex = L;
    end
  end

  % ------------------------- identify which bases each node goes with
      
  N = length(Node);
  d = Node(1).basehalfwidth;

  for n=1:N                                      % loop through nodes
    switch Node(n).type,
    case {'Initial','Basepair','Motif','Alternative','Junction'}
      iMin = max(min(Sequence(k).Constraint(n).LeftIndex)-d,1);
      iMax = min(max(Sequence(k).Constraint(n).LeftIndex)+d,L);
      jMin = max(min(Sequence(k).Constraint(n).RightIndex)-d,1);
      jMax = min(max(Sequence(k).Constraint(n).RightIndex)+d,L);

    case 'Hairpin'
      iMin = max(min(Sequence(k).Constraint(n).MiddleIndex)-d,1);
      iMax = min(min(Sequence(k).Constraint(n).MiddleIndex)+d,L);
      jMin = max(max(Sequence(k).Constraint(n).MiddleIndex)-d,1);
      jMax = min(max(Sequence(k).Constraint(n).MiddleIndex)+d,L);

    case {'JunctionCluster'}
      iMin = max(min(Sequence(k).Constraint(n).LeftIndex)-d,1);
      iMax = min(max(Sequence(k).Constraint(n).LeftIndex)+d,L);
      jMin = max(min(Sequence(k).Constraint(n).RightIndex)-d,1);
      jMax = min(max(Sequence(k).Constraint(n).RightIndex)+d,L);

      kMin = max(min(Sequence(k).Constraint(n).MiddleIndex)-d,1);
      kMax = min(max(Sequence(k).Constraint(n).MiddleIndex)+d,L);

      Sequence(k).Constraint(n).kMin = kMin;
      Sequence(k).Constraint(n).kMax = kMax;

    end

    Sequence(k).Constraint(n).iMin = iMin;
    Sequence(k).Constraint(n).iMax = iMax;
    Sequence(k).Constraint(n).jMin = jMin;
    Sequence(k).Constraint(n).jMax = jMax;

  end

end                                          % end loop over M sequences
    
%--------------------------------------------------------------------
% function V = pCount(V)
% 
% t=1;
% ok=0;
% while 1-ok
%     V(t)=V(t)+1;
%     ok=1;
%     if V(t)>4
%         V(t)=1;
%         t=t+1;
%         ok=0;
%     end
% end
