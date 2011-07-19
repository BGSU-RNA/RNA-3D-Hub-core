% pTheoreticalAlignment gives the header line and sequence that would
% result if each node generated exactly what it is designed to generate,
% as found in the 3D structure

function [Sequence] = pTheoreticalAlignment(Node,n)

Left{1}  = '';                  % header line
Left{2}  = '';                  % sequence line
Right{1} = '';
Right{2} = '';

while n <= length(Node),

%  fprintf('Node %d type %s\n', n, Node(n).type);

  switch Node(n).type
    case 'Initial' %-------------------------------------------------------
      Left{1}  = [Left{1} '['];
      Right{1} = [']' Right{1}];

      if length(Node(n).LeftLetter) > 1,
        for d = 2:length(Node(n).LeftLetter),
          Left{1} = [Left{1} '-'];
        end
        Left{2}  = [Left{2} Node(n).LeftLetter];
      elseif length(Node(n).LeftLetter) == 1,
        Left{2}  = [Left{2} Node(n).LeftLetter];
      else
        Left{2}  = [Left{2} '-'];
      end

      if length(Node(n).RightLetter) > 1,
        for d = 2:length(Node(n).RightLetter),
          Right{1} = ['-' Right{1}];
        end
        Right{2}  = [Node(n).RightLetter Right{2}];
      elseif length(Node(n).RightLetter) == 1,
        Right{2}  = [Node(n).RightLetter Right{2}];
      else
        Right{2}  = ['-' Right{2}];
      end

     

    case 'Basepair' % -----------------------------------------------------
      Left{1}  = [Left{1} '('];
      Right{1} = [')' Right{1}];

      if length(Node(n).LeftLetter) > 1,
        for d = 2:length(Node(n).LeftLetter),
          Left{1} = [Left{1} '-'];
        end
        Left{2}  = [Left{2} Node(n).LeftLetter];
      elseif length(Node(n).LeftLetter) == 1,
        Left{2}  = [Left{2} Node(n).LeftLetter];
      else
        Left{2}  = [Left{2} '-'];
      end

      if length(Node(n).RightLetter) > 1,
        for d = 2:length(Node(n).RightLetter),
          Right{1} = ['-' Right{1}];
        end
        Right{2}  = [Node(n).RightLetter Right{2}];
      elseif length(Node(n).RightLetter) == 1,
        Right{2}  = [Node(n).RightLetter Right{2}];
      else
        Right{2}  = ['-' Right{2}];
      end

    case 'Junction' %--------------------------------------------------------
      for c = 1:length(Node(n).nextnode),
        S       = pTheoreticalAlignment(Node,Node(n).nextnode(c));
        Left{1} = [Left{1} S{1}];
        Left{2} = [Left{2} S{2}];
      end
      n = Inf;


    case 'Cluster' % --------------------------------------------------------
      Left{1}  = [Left{1} '{'];
      a = '';
      for d = 1:length(Node(n).LeftLetter),
        a = [a '-'];
      end
      for d = 1:length(Node(n).Left),
        a(Node(n).Left(d)) = 'I';
      end
      Left{1} = [Left{1} a];

      Right{1} = ['}' Right{1}];
      a = '';
      for d = 1:length(Node(n).RightLetter),
        a = [a '-'];
      end
      for d = 1:length(Node(n).Right),
        a(Node(n).Right(d)) = 'I';
      end
      Right{1} = [a Right{1}];

      Left{2}  = [Left{2} '{' Node(n).LeftLetter];
      Right{2} = [Node(n).RightLetter '}' Right{2}];

    case 'Hairpin' % --------------------------------------------------------
      Left{1} = [Left{1} '<'];
      for d = 1:length(Node(n).LeftLetter),
        Left{1} = [Left{1} 'F'];
      end
      Right{1} = ['>' Right{1}];
      Left{2}  = [Left{2} '<' Node(n).LeftLetter '>'];
      n = Inf;
  end

  n = n + 1;

end

Sequence{1} = [Left{1} Right{1}];
Sequence{2} = [Left{2} Right{2}];

%fprintf('%s\n',Sequence{1});
%fprintf('%s\n',Sequence{2});
%fprintf('\n');