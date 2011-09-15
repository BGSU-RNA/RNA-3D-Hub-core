% pAlignSequence builds the alignment of a single sequence

function [alignment] = pAlignSequence(Node,t,maxinsert,n,s,i,j,x)

left     = '';
right    = '';
middle   = '';
AltNode=0;
N=0;

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
       Letters(pAlignSequence(Node,t,maxinsert,Node(n).nextnode(1),s,i,a-1,x),active)];
     right = [...
       Letters(pAlignSequence(Node,t,maxinsert,Node(n).nextnode(2),s,a,j,x),active) right];
   
    case 'Alternative'
        % see below
        
    case 'Basepair'
      if (s==1),
        left  = [left Char(active,maxinsert(n).left+1)];
        right = [Char(active,maxinsert(n).right+1) right];
      else
           if active
            left = [left Letters(x(i:(a-1)),active) Char(active,maxinsert(n).left-(a-i-1))];
            right = [Char(active,maxinsert(n).right-(j-b-1)) Letters(x((b+1):j),active) right];
           else
%             left = [left Letters(x(i:i),active) Char(active,maxinsert(n).left-(a-i-1))];
%             right = [Char(active,maxinsert(n).right-(j-b-1)) Letters(x(j:j),active) right];
            left = [left Letters('A',active) Char(active,maxinsert(n).left)];
            right = [Char(active,maxinsert(n).right) Letters('A',active) right];
           end
      end
      
    case 'JunctionCluster' 
     m    = double(t(n).m);
     c    = double(t(n).c);
     l    = double(t(n).l);
     r    = double(t(n).r);
     cm   = double(t(n).cm);
     left = ClusterLMR(left,Node(n).Left,maxinsert(n).left,n,l,x,i,active,1);
     left = [left ...
       Letters(pAlignSequence(Node,t,maxinsert,Node(n).nextnode(1),s,i,m,x),active)];
     middle=ClusterLMR(middle,Node(n).Middle,maxinsert(n).middle,n,cm,x,m+1,active,1);
     right=ClusterLMR(right,Node(n).Right,maxinsert(n).right,n,r,x,j,active,-1);
     right = [...
             Letters(pAlignSequence(Node,t,maxinsert,Node(n).nextnode(2),s,c,j,x),active) right];       
     
    case 'Cluster'
     l = double(t(n).l);
     r = double(t(n).r);
     if (s==1)                   % deleted
       left = [left Char(active,length(Node(n).Left(1,:)) + sum(maxinsert(n).left))];
       right= [Char(active,length(Node(n).Right(1,:))+sum(maxinsert(n).right)) right];
     else
      if t(n).active==1
          left = ClusterLMR(left ,Node(n).Left ,maxinsert(n).left ,n,l,x,i,active, 1);
          right= ClusterLMR(right,Node(n).Right,maxinsert(n).right,n,r,x,j,active,-1);
      else
          L1=length(Node(n).Left(1,:));
          R1=length(Node(n).Right(1,:));
          
          left =[left Char(active,L1) Char(active,maxinsert(n).left)];
          right= [Char(active,maxinsert(n).right) Char(active,R1) right];
      end
     end

  end
%   n
%   active
    if AltNode>0 & AltNode<length(N)
        if Node(n).nextnode(1)>N(AltNode+1);
            n = N(AltNode+1);
            AltNode=AltNode+1;
        else
            [n,N,AltNode]=Switch(n,N,AltNode,Node);
        end
    else 
        [n,N,AltNode]=Switch(n,N,AltNode,Node);
    end
%     AltNode
%     
%     disp('--------------')
end

alignment = [left middle right];
%----------------------------------------------------------------
function [position]=ClusterLMR(position,NodeP,maxinsertP,n,l,x,i,active,p)
    if p==1
       A = 1:length(NodeP(1,:))-1; 
       q=0;
    else  
       A = length(NodeP(1,:)):-1:2;
       q=-1;
    end
      A = 1:length(NodeP(1,:))-1; 
      bb=0;
    for k=1:length(A)
        a=A(k);
        if q==-1
            position = [Letters(x(i+1-a-bb),active) position];
            for b=1:maxinsertP(a),
              if b < NodeP(l,a) - NodeP(l,a+1),
                  bb=bb+1;          
                position = [Letters(x(i+1-a-bb),active) position];
              else
                position = [Char(active) position];
              end
            end
        else  
            position = [position Letters(x(i-1+NodeP(l,a)),active)];
            for b=1:maxinsertP(a),
                
              if b < NodeP(l,a+1) - NodeP(l,a),
                  
                position = [position Letters(x(i-1+NodeP(l,a)+b),active)];
              else
                position = [position Char(active)];
              end
            end
        end
    end
    if p==1
       position = [position Letters(x(i-1+NodeP(l,length(NodeP(1,:)))),active)];
    else
       position = [Letters(x(i+1-NodeP(l,1)),active) position];
    end

%--------------------------------------------------------------------
    function position=Order(x,y,p)
        if p==1;
            position=[x y];
        else
            position=[y x];
        end
 % -----------------------------------------  
function [n,N,AltNode]=Switch(n,N,AltNode,Node)
    switch Node(n).type,
    case {'Initial', 'Basepair', 'Cluster'}
      n = Node(n).nextnode(1);
    case 'Alternative'
        N=Node(n).nextnode;
        AltNode=1;
        n=Node(n).nextnode(1); 
    case {'Hairpin', 'Junction','JunctionCluster'}
      n = Inf;
    end
 % -----------------------------------------    
function letters=Letters(y,active)
    if active==1
        letters=y;
    else
        letters=Char(0,length(y));
    end     
% -------------------------- -'s or +'s of a specified length     
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
