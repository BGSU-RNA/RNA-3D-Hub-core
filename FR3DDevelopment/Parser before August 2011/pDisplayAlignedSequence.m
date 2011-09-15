% pDisplayAlignedSequence displays one sequence in a multiple alignment

function [alignment] = pDisplayAlignedSequence(Node,t,maxinsert,n,s,i,j,x)

left  = '';
right = '';

while (n <= length(Node)),                          % go through nodes
  i = double(t(n).i);
  j = double(t(n).j);
  a = double(t(n).a);
  b = double(t(n).b);
  s = t(n).state;

  switch Node(n).type,
    case 'Hairpin'
      left = [left '*'];
      right = ['*' right];
      left = [left x(i:j)];
      for a=1:(maxinsert(n).hairpin-(j-i+1)),
        left = [left '-'];
      end
    case 'Initial'
      left = [left x(i:(a-1))];
      right = [x((b+1):j) right];
      for a=1:(maxinsert(n).left-(a-i)),
        left = [left '-'];
      end
      for a=1:(maxinsert(n).right-(j-b)),
        right = ['-' right];
      end
    case 'Junction'
     left = [left ...
       pDisplayAlignedSequence(Node,t,maxinsert,Node(n).nextnode(1),s,i,a-1,x)];
     right = [...
       pDisplayAlignedSequence(Node,t,maxinsert,Node(n).nextnode(2),s,a,j,x) ...
       right];
    case 'Generic'
      if (s==1),
        for a=1:maxinsert(n).left+1,
          left = [left '-'];
        end
        for a=1:maxinsert(n).right+1,
          right = ['-' right];
        end
      else
        left = [left x(i:(a-1))];
        right = [x((b+1):j) right];
        for a=1:(maxinsert(n).left-(a-i-1)),
          left = [left '-'];
        end
        for a=1:(maxinsert(n).right-(j-b-1)),
          right = ['-' right];
        end
      end
    case 'Motif'
  l = double(t(n).l);
  r = double(t(n).r);
     if (s==1),                                 % deleted
       for a=1:(length(Node(n).Left(1,:)) + sum(maxinsert(n).left)),
         left = [left '-'];
       end
       for a=1:(length(Node(n).Right(1,:)) + sum(maxinsert(n).right)),
         right = ['-' right];
       end
     else
      for a = 1:length(Node(n).Left(1,:))-1,
        left = [left x(i-1+Node(n).Left(l,a))];
        for b=1:maxinsert(n).left(a),
          if b < Node(n).Left(l,a+1) - Node(n).Left(l,a),
            left = [left x(i-1+Node(n).Left(l,a)+b)];
          else
            left = [left '-'];
          end
        end
      end
      left = [left x(i-1+Node(n).Left(l,length(Node(n).Left(1,:))))];

      for a = length(Node(n).Right(1,:)):-1:2,
        right = [x(j+1-Node(n).Right(r,a)) right];
        for b=1:maxinsert(n).right(a-1),
          if b < Node(n).Right(r,a-1) - Node(n).Right(r,a),
            right = [x(j+1-Node(n).Right(r,a)-b) right];
          else
            right = ['-' right];
          end
        end
      end
      right = [x(j+1-Node(n).Right(r,1)) right];
     end
  end

  switch Node(n).type,
    case {'Initial', 'Generic', 'Motif'}
      n = Node(n).nextnode(1);
    case {'Hairpin', 'Junction'}
      n = Inf;
  end

end

alignment = [left right];
