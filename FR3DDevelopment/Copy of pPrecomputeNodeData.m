% pPrecomputeNodeData(Node) computes logarithms and such that will be
% needed when parsing according to the model in Node

function [Node] = pPrecomputeNodeData(Node)

% ---------------------------------- Precompute logarithms, etc. for Node

pairs = {'AA', 'CA', 'GA', 'UA', 'AC', 'CC', 'GC', 'UC', 'AG', 'CG', 'GG', 'UG', 'AU', 'CU', 'GU', 'UU', 'DEL'};

k = 10;                     % maximum number of insertions to calculate for

warning off

for n=1:N,
  switch Node(n).type,
    case 'Initial'
      Node(n).numstates = 1;
    case 'Basepair'
      Node(n).numstates = 2;    % deleted or active
    case 'Junction'
      Node(n).numstates = 1;
    case 'Hairpin'
     Node(n).numstates = 1;
    case 'Motif'
      Node(n).numstates = 2;    % deleted or active
    case 'Alternative'
      Node(n).numstates = 1;
  end

  switch Node(n).type,
    case {'Basepair','Hairpin','Junction','Motif'}
      Node(n).lP         = log(Node(n).P);           
      Node(n).lPIns      = log(Node(n).PIns);           
  end

  switch Node(n).type,
    case {'Motif'}
      Node(n).LLIP       = log(Node(n).LIP);           
      Node(n).LRIP       = log(Node(n).RIP);           
  end

  switch Node(n).type,
    case {'Initial','Basepair'}
      for s=1:length(Node(n).lpar(:,1)),
        l = Node(n).lpar(s);
        Node(n).lip(s,:) = exp(-l)*((ones(1,k+1)*l).^(0:k))./[1 cumprod(1:k)];
        r = Node(n).rpar(s);
        Node(n).rip(s,:) = exp(-r)*((ones(1,k+1)*r).^(0:k))./[1 cumprod(1:k)];
      end

      Node(n).loglip = log(Node(n).lip);
      Node(n).logrip = log(Node(n).rip);
  end
end

warning on

fprintf('Did initial internal calculations\n');
