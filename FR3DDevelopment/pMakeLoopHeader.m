% pMakeLoopHeader constructs a 3-line header for a multiple alignment

function [header] = pMakeLoopHeader(Node,maxinsert,n)

left1  = '';
left2  = '';
left3  = '';

right1 = '';
right2 = '';
right3 = '';
AltNode=0;
N=1;
while (n <= length(Node)),                            % go through nodes
%   active=t(n).active;
  switch Node(n).type,
    case 'Hairpin'
      left1  = [left1 '*H'];
      left2  = [left2 '*' num2str(mod(n,10))];
      left3  = [left3 '*' num2str(fix(n/10))];

      for a=1:maxinsert(n).hairpin-1,
        left1 = [left1 '-'];
        left2 = [left2 num2str(mod(n,10))];
        left3 = [left3 num2str(fix(n/10))];
      end

      right1 = ['*' right1];
      right2 = ['*' right2];
      right3 = ['*' right3];
    case 'Alternative'
        
    
    case 'Initial'
      for a=1:maxinsert(n).left,
        left1 = [left1 '-'];
        left2 = [left2 num2str(mod(n,10))];
        left3 = [left3 num2str(fix(n/10))];
      end
      for a=1:maxinsert(n).right,
        right1 = ['-' right1];
        right2 = [num2str(mod(n,10)) right2];
        right3 = [num2str(fix(n/10)) right3];
      end

    case 'Junction'
      t = pMakeLoopHeader(Node,maxinsert,Node(n).nextnode(1));
      left1  = [left1 t{1}];
      left2  = [left2 t{2}];
      left3  = [left3 t{3}];

      t = pMakeLoopHeader(Node,maxinsert,Node(n).nextnode(2));
      right1 = [ t{1} right1];
      right2 = [ t{2} right2];
      right3 = [ t{3} right3];
      
    case 'Basepair'
      left1  = [left1 '('];
      left2  = [left2 num2str(mod(n,10))];
      left3  = [left3 num2str(fix(n/10))];

      for a=1:maxinsert(n).left,
        left1 = [left1 '-'];
        left2 = [left2 num2str(mod(n,10))];
        left3 = [left3 num2str(fix(n/10))];
      end

      right1 = [')' right1];
      right2 = [num2str(mod(n,10)) right2];
      right3 = [num2str(fix(n/10)) right3];

      for a=1:maxinsert(n).right,
        right1 = ['-' right1];
        right2 = [num2str(mod(n,10)) right2];
        right3 = [num2str(fix(n/10)) right3];
      end
      
    case 'JunctionCluster'
       % start Left -------------------------------
      for a = 1:length(Node(n).Left(1,:))-1,
        left1 = [left1 'J'];
        for b=1:maxinsert(n).left(a),
          left1 = [left1 'B'];
        end
      end
      left1 = [left1 'J'];
      for a=1:(length(Node(n).Left(1,:)) + sum(maxinsert(n).left)),
        left2  = [left2 num2str(mod(n,10))];
        left3  = [left3 num2str(fix(n/10))];
      end  
      % end Left -------------------------------
      %  1st helix------------------------------
      t = pMakeLoopHeader(Node,maxinsert,Node(n).nextnode(1));
      left1  = [left1 t{1}];
      left2  = [left2 t{2}];
      left3  = [left3 t{3}];
      % start Middle -------------------------------
      for a = 1:length(Node(n).Middle(1,:))-1,
        left1 = [left1 'J'];
        for b=1:maxinsert(n).middle(a),
          left1 = [left1 'B'];
        end
      end
      left1 = [left1 'J'];
      for a=1:(length(Node(n).Middle(1,:)) + sum(maxinsert(n).middle)),
        left2  = [left2 num2str(mod(n,10))];
        left3  = [left3 num2str(fix(n/10))];
      end     
      % end Middle -------------------------------
      % start Right -------------------------------
      for a = 1:length(Node(n).Right(1,:))-1,
        right1 = ['J' right1];
        for b=1:maxinsert(n).right(a),
          right1 = ['B' right1];
        end
      end
      right1 = ['J' right1];

      for a=1:(length(Node(n).Right(1,:)) + sum(maxinsert(n).right)),
        right2 = [num2str(mod(n,10)) right2];
        right3 = [num2str(fix(n/10)) right3];
      end
      % end Right -------------------------------
      %  2nd helix------------------------------
      t = pMakeLoopHeader(Node,maxinsert,Node(n).nextnode(2));
      right1 = [ t{1} right1];
      right2 = [ t{2} right2];
      right3 = [ t{3} right3];    
      
      
    case 'Cluster'
      for a = 1:length(Node(n).Left(1,:))-1,  % loop through left letters
        left1 = [left1 'C'];
[n a Node(n).Left(1,:)]
maxinsert(n).left
        for b=1:maxinsert(n).left(a),
          left1 = [left1 'B'];
        end
      end
      left1 = [left1 'C'];
      for a = 1:length(Node(n).Right(1,:))-1,
        right1 = ['C' right1];
        for b=1:maxinsert(n).right(a),
          right1 = ['B' right1];
        end
      end
      right1 = ['C' right1];

      for a=1:(length(Node(n).Left(1,:)) + sum(maxinsert(n).left)),
        left2  = [left2 num2str(mod(n,10))];
        left3  = [left3 num2str(fix(n/10))];
      end
      for a=1:(length(Node(n).Right(1,:)) + sum(maxinsert(n).right)),
        right2 = [num2str(mod(n,10)) right2];
        right3 = [num2str(fix(n/10)) right3];
      end
  end

  switch Node(n).type,
    case 'Alternative'
        N=Node(n).nextnode;
        AltNode=1;
        n=Node(n).nextnode(1);
    case {'Initial', 'Basepair', 'Cluster'}
         n = Node(n).nextnode(1);
       if AltNode>0 & AltNode<length(N)
           if n>N(AltNode+1);
               n = N(AltNode+1);
               AltNode=AltNode+1;
           end
       end
    case {'Hairpin', 'Junction','JunctionCluster'}
      n = Inf;
  end

%fprintf('pMakeLoopHeader Node %d\n', n);
%fprintf('%s %s\n', left3, right3);
%fprintf('%s %s\n', left2, right2);
%fprintf('%s %s\n', left1, right1);

end

header{1} = [left1 right1];
header{2} = [left2 right2];
header{3} = [left3 right3];

