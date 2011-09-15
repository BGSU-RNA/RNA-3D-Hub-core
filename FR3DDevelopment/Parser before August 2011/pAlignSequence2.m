% pAlignSequence2 builds the alignment of a single sequence

function [alignment] = pAlignSequence2(Node,t,maxinsert,n,s,i,j,x)
% x
left        = '';
right      = '';
middle   = '';
% 
%   for n=1:length(Node)
%       t(n).state
%   end
  AltNode=0;
  N=0;
%   Char(0,1)
while (n <= length(Node)),                          % go through nodes

  i = double(t(n).i);
  j = double(t(n).j);
  a = double(t(n).a);
  b = double(t(n).b);
  s = t(n).state;
  active=t(n).active;

%   [n active]
  switch Node(n).type,
    case 'Hairpin'
      left = [left '*' Letters(x(i:j),active) Char(active,maxinsert(n).hairpin-(j-i+1))];
      right = ['*' right];
    case 'Initial'
      left = [left Letters(x(i:(a-1)),active) Char(active,maxinsert(n).left-(a-i))];
      right = [Char(active,maxinsert(n).right-(j-b)) x((b+1):j) right];
    case 'Junction'
     left = [left ...
       Letters(pAlignSequence22(Node,t,maxinsert,Node(n).nextnode(1),s,i,a-1,x),active)];
     right = [...
       Letters(pAlignSequence22(Node,t,maxinsert,Node(n).nextnode(2),s,a,j,x),active) right];
    case 'Alternative'
        % see below
    case 'Basepair'
%         if n==36
%             active
%             Char(active,maxinsert(n).left+1)
%             [Letters(x(i:(a-1)),active) Char(active,maxinsert(n).left-(a-i-1))]
%             [Char(active,maxinsert(n).right-(j-b-1)) Letters(x((b+1):j),active)]
%         end
        
      if (s==1),
        left  = [left Char(active,maxinsert(n).left+1)];
        right = [Char(active,maxinsert(n).right+1) right];
      else
        left = [left Letters(x(i:(a-1)),active) Char(active,maxinsert(n).left-(a-i-1))];
        right = [Char(active,maxinsert(n).right-(j-b-1)) Letters(x((b+1):j),active) right];
      end
    case 'JunctionMotif' 
     m = double(t(n).m);
     c = double(t(n).c);
     l  = double(t(n).l);
     r  = double(t(n).r);
     cm = double(t(n).cm);
     if (s==2),                                 % deleted

%        left = [left Char(active,length(Node(n).Left(1,:)) + sum(maxinsert(n).left))];
%        right= [Char(active,length(Node(n).Right(1,:))+sum(maxinsert(n).right)) right];
     else
      for a = 1:length(Node(n).Left(1,:))-1,
        left = [left Letters(x(i-1+Node(n).Left(l,a)),active)];
        for b=1:maxinsert(n).left(a),
          if b < Node(n).Left(l,a+1) - Node(n).Left(l,a),
            left = [left Letters(x(i-1+Node(n).Left(l,a)+b),active)];
          else
            left = [left Char(active,1)];
          end % end if 
        end % end b
      end % end a
      left = [left Letters(x(i-1+Node(n).Left(l,length(Node(n).Left(1,:)))),active)];

     left = [left ...
       Letters(pAlignSequence2(Node,t,maxinsert,Node(n).nextnode(1),s,i,m,x),active)];
   % end left part
   % begining Middle
      for a = 1:length(Node(n).Middle(1,:))-1,
        middle = [middle Letters(x(m+Node(n).Middle(cm,a)),active)];
        for b=1:maxinsert(n).middle(a),
          if b < Node(n).Middle(cm,a+1) - Node(n).Middle(cm,a),
            middle = [middle Letters(x(m+Node(n).Middle(cm,a)+b),active)];
          else
            middle = [middle Char(active,1)];
          end % end if 
        end % end b
      end % end a
      middle = [middle Letters(x(m+Node(n).Middle(cm,length(Node(n).Middle(1,:)))),active)];
% end Middle

      for a = length(Node(n).Right(1,:)):-1:2,
        right = [Letters(x(j+1-Node(n).Right(r,a)),active) right];
        for b=1:maxinsert(n).right(a-1),
          if b < Node(n).Right(r,a-1) - Node(n).Right(r,a),
            right = [Letters(x(j+1-Node(n).Right(r,a)-b),active) right];
          else
            right = [Char(active,1) right];
          end % end if 
        end % end b
      end % end a
      right = [Letters(x(j+1-Node(n).Right(r,1)),active) right];
      right = [...
             Letters(pAlignSequence2(Node,t,maxinsert,Node(n).nextnode(2),s,c,j,x),active) right];       
     end

    case 'Motif'
     l = double(t(n).l);
     r = double(t(n).r);
     if length(l)==0 & active==0
         l=1;
         r=1;
         i=1;
         j=length(x);
     end
     if (s==1)                   % deleted
       left = [left Char(active,length(Node(n).Left(1,:)) + sum(maxinsert(n).left))];
       right= [Char(active,length(Node(n).Right(1,:))+sum(maxinsert(n).right)) right];
     else
      for a = 1:length(Node(n).Left(1,:))-1,
        left = [left Letters(x(i-1+Node(n).Left(l,a)),active)];
        [n a+1 length(Node(n).Left(1,:))-1]
         for b=1:maxinsert(n).left(a+1), 
          if b < Node(n).Left(l,a+1) - Node(n).Left(l,a),
            left = [left Letters(x(i-1+Node(n).Left(l,a)+b),active)];
          else
            left = [left Char(active)];
          end
        end
      end
      left = [left Letters(x(i-1+Node(n).Left(l,length(Node(n).Left(1,:)))),active)];
      for a = length(Node(n).Right(1,:)):-1:2,
        right = [Letters(x(j+1-Node(n).Right(r,a)),active) right];
        for b=1:maxinsert(n).right(a-1),
          if b < Node(n).Right(r,a-1) - Node(n).Right(r,a),
            right = [Letters(x(j+1-Node(n).Right(r,a)-b),active) right];
          else
            right = [Char(active) right];
          end
        end
      end
      right = [Letters(x(j+1-Node(n).Right(r,1)),active) right];
     end

  end
    if AltNode>0 & AltNode<length(N)
    %             [n N(AltNode+1) 111]
        if Node(n).nextnode(1)>N(AltNode+1);
            n = N(AltNode+1);
            AltNode=AltNode+1;
        else
           switch Node(n).type,
                case {'Initial', 'Basepair', 'Motif'}
                  n = Node(n).nextnode(1);
                case 'Alternative'
                    N=Node(n).nextnode;
                    AltNode=1;
                    n=Node(n).nextnode(1); 
                case {'Hairpin', 'Junction','JunctionMotif'}
                  n = Inf;
           end       
        end
    else 
      switch Node(n).type,

        case {'Initial', 'Basepair', 'Motif'}
          n = Node(n).nextnode(1);
        case 'Alternative'
            N=Node(n).nextnode;
            AltNode=1;
            n=Node(n).nextnode(1); 
        case {'Hairpin', 'Junction','JunctionMotif'}
          n = Inf;
      end
    end



end

alignment = [left middle right];

% -------------------------- Dashes of a specified length


     
function letters=Letters(y,active)
    if active==1
        letters=y;
    else
        letters=Char(0,length(y));
    end     
     
function [ch]=Char(k,n)
    if nargin==1
        n=1;
    end
    switch k
        case 0
            ch=char('+'*ones(1,n));
        case 1
            ch=char('-'*ones(1,n));
    end



function [dashes] = Dash(n)

    dashes = char('-'*ones(1,n));

function [plusses] = Plus(n)

    dashes = char('+'*ones(1,n));

