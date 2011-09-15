% pAlignSequence3 builds the alignment of a single sequence

function [alignment] = pAlignSequence3(Node,t,maxinsert,n,s,i,j,x)

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
      left = [left '*' x(i:j) Dash(maxinsert(n).hairpin-(j-i+1))];
      right = ['*' right];
    case 'Initial'
      left = [left x(i:(a-1)) Dash(maxinsert(n).left-(a-i))];
      right = [Dash(maxinsert(n).right-(j-b)) x((b+1):j) right];
    case 'Junction'
     left = [left ...
       pAlignSequence3(Node,t,maxinsert,Node(n).nextnode(1),s,i,a-1,x)];
     right = [...
       pAlignSequence3(Node,t,maxinsert,Node(n).nextnode(2),s,a,j,x) right];
    case 'Basepair'
      if (s==1),
        left  = [left Dash(maxinsert(n).left+1)];
        right = [Dash(maxinsert(n).right+1) right];
      else
        left = [left x(i:(a-1)) Dash(maxinsert(n).left-(a-i-1))];
        right = [Dash(maxinsert(n).right-(j-b-1)) x((b+1):j) right];
      end
    case 'Motif'
     l = double(t(n).l);
     r = double(t(n).r);
     if (s==1),                                 % deleted
       left = [left Dash(length(Node(n).Left(1,:)) + sum(maxinsert(n).left))];
       right= [Dash(length(Node(n).Right(1,:))+sum(maxinsert(n).right)) right];
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
    case {'Initial', 'Basepair', 'Motif'}
      n = Node(n).nextnode(1);
    case {'Hairpin', 'Junction'}
      n = Inf;
  end

end

alignment = [left right];

% -------------------------- Dashes of a specified length

function [dashes] = Dash(n)

dashes = char('-'*ones(1,n));

