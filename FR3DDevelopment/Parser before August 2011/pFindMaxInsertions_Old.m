% pFindMaxInsertions(Node,t,maxinsert,L) determines the largest number of
% insertions, among all sequences, for each possible location of insertions
% in Node

function [maxinsert] = pFindMaxInsertions(Node,t,maxinsert,L)

N = length(Node);                         % number of nodes

for n=1:N,
  i = double(t(n).i);
  j = double(t(n).j);
  a = double(t(n).a);
  b = double(t(n).b);
  s = t(n).state;
  active=t(n).active;
%   [n active]
  switch Node(n).type,
    case 'Hairpin'
      if active
        maxinsert(n).hairpin = max(maxinsert(n).hairpin,j-i+1);
      else
        maxinsert(n).hairpin = max(maxinsert(n).hairpin,1);
      end
    case 'Initial'  ;
      if active
        maxinsert(n).left  = max(maxinsert(n).left,  a-i);
        maxinsert(n).right = max(maxinsert(n).right, j-b);
      else
        maxinsert(n).left  = max(maxinsert(n).left,0);
        maxinsert(n).right = max(maxinsert(n).right,0);
      end
    case 'Basepair'
      if active
        maxinsert(n).left  = max(maxinsert(n).left, a-i-1);
        maxinsert(n).right = max(maxinsert(n).right,j-b-1);
      else
        maxinsert(n).left  = max(maxinsert(n).left,0);
        maxinsert(n).right = max(maxinsert(n).right,0);  
      end
    case {'Motif'}
      if active
          if (s==2),                         % motif is not deleted
            l = double(t(n).l);              % left  insert pattern used in motif
            r = double(t(n).r);              % right insert pattern used in motif
            d = diff(Node(n).Left(l,:));          % differences used

            maxinsert(n).left  = max(maxinsert(n).left, d-1);
            d = -diff(Node(n).Right(r,:));          % differences used
            maxinsert(n).right = max(maxinsert(n).right, d-1);
          else                               % motif is deleted
            maxinsert(n).left  = max(maxinsert(n).left, 0);
            maxinsert(n).right = max(maxinsert(n).right,0);
          end
      else
          
          maxinsert(n).left  = max(maxinsert(n).left, 0);
          maxinsert(n).right = max(maxinsert(n).right,0);
      end
    case {'JunctionMotif'}
      if active
        l = double(t(n).l);              % left  insert pattern used in motif
        r = double(t(n).r);              % right insert pattern used in motif
        c = double(t(n).cm);             % middle insert pattern used in motif
        d = diff(Node(n).Left(l,:));          % differences used
        maxinsert(n).left  = max(maxinsert(n).left, d-1);
        d = diff(Node(n).Middle(c,:));          % differences used
        maxinsert(n).middle  = max(maxinsert(n).middle, d-1);        
        d = -diff(Node(n).Right(r,:));         % differences used
        maxinsert(n).right = max(maxinsert(n).right, d-1);  
      else
        maxinsert(n).left   = max(maxinsert(n).left, 0);
        maxinsert(n).middle = max(maxinsert(n).middle,0);       
        maxinsert(n).right  = max(maxinsert(n).right,0);   
      end
  end
end

