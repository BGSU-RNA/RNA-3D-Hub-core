% pWriteJavaNodeFile(File,Node) writes Java code describing each JAR3D node

% NumChar tells it how large to make letter and pair distributions, what to
% add for * hairpins, for example

function [] = pWriteJavaNodeFile(File,Node,NumChar,Filename)

if nargin < 3,
  NumChar = 4;
end

if nargin < 4,
  fid = fopen(['Models/' File.Filename '.txt'],'w');
else
  fid = fopen(['Models/' Filename], 'w');
end

if NumChar == 4,
  Text = 'Character Definition | A,C,G,U // Define characters here';
elseif NumChar == 5,
  Text = 'Character Definition | A,C,G,U,* // Define characters here';
elseif NumChar == 6,
  Text = 'Character Definition | A,C,G,U,N // Define characters here';
end
fprintf(fid,'%s\n',Text);

for n=1:length(Node),
  Text = '';

  switch Node(n).type
    case 'Initial' % -----------------------------------------------------

      Text = [Text sprintf('InitialNode  | ')];

      Text = [Text sprintf('Left Length Distribution ')];
      Text = [Text subWrite(fid, Node(n).leftLengthDist)];

      Text = [Text sprintf(' | Left Letter Distribution ')];
      Text = [Text subWrite(fid, Node(n).leftLetterDist, NumChar)];

      Text = [Text sprintf(' | Right Length Distribution ')];
      Text = [Text subWrite(fid, Node(n).rightLengthDist)];

      Text = [Text sprintf(' | Right Letter Distribution ')];
      Text = [Text subWrite(fid, Node(n).rightLetterDist, NumChar)];

      LI = Node(n).LeftIndex;
      RI = Node(n).RightIndex;

      Text = [Text sprintf(' | Left Index [%d]', LI(1))];
      Text = [Text sprintf(' | Right Index [%d]', RI(end))];

      Text = [Text Node(n).Comment sprintf('\n')];

    case 'Junction' % -----------------------------------------------------
      Text = [Text sprintf('//\nJunctionNode | Branches ')];
      Text = [Text sprintf('[%d] ',length(Node(n).nextnode))];

      LI = Node(n).LeftIndex;
      RI = Node(n).RightIndex;

      Text = [Text sprintf(' | Left Index [%d]', LI(1))];
      Text = [Text sprintf(' | Right Index [%d] ', RI(end))];

      Text = [Text Node(n).Comment sprintf('\n')];

    case 'Basepair' % ------------------------------------------------------

      Text = [Text sprintf('BasepairNode | Deletion Probability ')];
      Text = [Text sprintf('[%0.9f] | ', Node(n).Delete)];

      Text = [Text sprintf('Pair Probability ')];
      Text = [Text subWriteArray(fid, Node(n).SubsProb, NumChar)];

      Text = [Text sprintf(' | Left Length Distribution ')];
      Text = [Text subWrite(fid, Node(n).leftLengthDist)];

      Text = [Text sprintf(' | Left Letter Distribution ')];
      Text = [Text subWrite(fid, Node(n).leftLetterDist, NumChar)];

      Text = [Text sprintf(' | Right Length Distribution ')];
      Text = [Text subWrite(fid, Node(n).rightLengthDist)];

      Text = [Text sprintf(' | Right Letter Distribution ')];
      Text = [Text subWrite(fid, Node(n).rightLetterDist, NumChar)];

      LI = Node(n).LeftIndex;
      RI = Node(n).RightIndex;

      Text = [Text sprintf(' | Left Index [%d]', LI(1))];
      Text = [Text sprintf(' | Right Index [%d]', RI(end))];

      Text = [Text Node(n).Comment sprintf('\n')];

    case 'Cluster' % -----------------------------------------------------

      Text = [Text sprintf('ClusterNode  | Deletion Probability ')];
      Text = [Text sprintf('[%0.9f] | ', Node(n).Delete)];
      Text = [Text sprintf('Num Left Bases [%d] | ', length(Node(n).Left))];
      Text = [Text sprintf('Num Right Bases [%d]', length(Node(n).Right))];

      LI = Node(n).LeftIndex;                      % indices on the left
      RI = Node(n).RightIndex;                     % indices on the right

      Text = [Text sprintf(' | Left Index [%d]', LI(1))];
      Text = [Text sprintf(' | Right Index [%d]', RI(end))];
      Text = [Text sprintf(' | Norm Constant [%0.15f]', Node(n).NormCons)];
      Text = [Text Node(n).Comment sprintf('\n')];

      % add interactions between bases ------------------------------------

      Indices = [Node(n).LeftIndex(Node(n).Left) ...
                 Node(n).RightIndex(Node(n).Right)];

      for i = 1:length(Node(n).IBases(:,1)),
        Text = [Text sprintf('Interaction  | Interacting Bases ')];
        Text = [Text sprintf('[%d,%d] ', ...
                     Node(n).IBases(i,1), Node(n).IBases(i,2))];

        Text = [Text sprintf('| Pair Probability ')];

        Text = [Text subWriteArray(fid, Node(n).SubsProb(:,:,i) , NumChar)];

        Text = [Text Node(n).InteractionComment{i} sprintf('\n')];
      end

      % add insertions between bases ---------------------------------------

      if ~isempty(Node(n).Insertion),
       for i = 1:length(Node(n).Insertion),
         p = Node(n).Insertion(i).Position;
         Text = [Text sprintf('Insertion    | Location [%d] ',p)];

         Text = [Text sprintf('| Length Distribution ')];
         Text = [Text subWrite(fid, Node(n).Insertion(i).LengthDist)];

         Text = [Text sprintf(' | Letter Distribution ')];
         Text = [Text subWrite(fid, Node(n).Insertion(i).LetterDist, NumChar)];

         Text = [Text Node(n).InsertionComment{i} sprintf('\n')];
       end
      end

    case 'Hairpin' % ------------------------------------------------------

      Text = [Text sprintf('HairpinNode  | ')];
      Text = [Text sprintf('Num Bases [%d]', length(Node(n).MiddleIndex))];
      
      if numel(Node(n).MiddleIndex) > 0 && Node(n).MiddleIndex(1) ~= -1,
        LI = Node(n).MiddleIndex(1);                      % use the first
        RI = Node(n).MiddleIndex(end);                    % use the last
      else
        LI = Node(n).LeftIndex(1);                      % use the first
        RI = Node(n).RightIndex(end);                    % use the last  
        fprintf('pWriteJavaNodeFile: Warning, Node(%d).MiddleIndex is empty\n',n);
      end

      Text = [Text sprintf(' | Left Index [%d]', LI)];
      Text = [Text sprintf(' | Right Index [%d]', RI)];
      Text = [Text sprintf(' | Norm Constant [%0.15f]', Node(n).NormCons)];
      Text = [Text Node(n).Comment sprintf('\n')];

      % add interactions between bases ------------------------------------

      if isfield(Node(n),'IBases'),
       if ~isempty(Node(n).IBases),
        for i = 1:length(Node(n).IBases(:,1)),
          Text = [Text sprintf('Interaction  | Interacting Bases ')];
          Text = [Text sprintf('[%d,%d] ', ...
                       Node(n).IBases(i,1), Node(n).IBases(i,2))];

          Text = [Text sprintf('| Pair Probability ')];

          Text = [Text subWriteArray(fid, Node(n).SubsProb(:,:,i), NumChar)];

          Text = [Text Node(n).InteractionComment{i} sprintf('\n')];
        end
       end
      end

      % add insertions between bases ---------------------------------------

      if ~isempty(Node(n).Insertion),
       for i = 1:length(Node(n).Insertion),
         p = Node(n).Insertion(i).Position;
         Text = [Text sprintf('Insertion    | Location [%d] ',p)];

         Text = [Text sprintf('| Length Distribution ')];
         Text = [Text subWrite(fid, Node(n).Insertion(i).LengthDist)];

         Text = [Text sprintf(' | Letter Distribution ')];
         Text = [Text subWrite(fid, Node(n).Insertion(i).LetterDist, NumChar)];

         Text = [Text Node(n).InsertionComment{i} sprintf('\n')];
       end
      end

  end

  fprintf(fid,'%s',Text);

end

fclose(fid);

% ---------------------------------------------------------------------------

function [newText] = subWrite(fid, a, L)

  newText = '';
  newText = [newText sprintf('[')];
  a2(1:length(a))=a;
  if nargin > 2,
    if length(a) < L,
      z = zeros(1,L-length(a));
      a2(length(a)+1:length(a)+length(z)) = z;
    end
  end
  a = a2;
  for i = 1:length(a)-1,
    newText = [newText sprintf('%0.9f,', a(i))];
  end
  newText = [newText sprintf('%0.9f]', a(end))];

% ---------------------------------------------------------------------------

function [newText] = subWriteArray(fid, a, NumChar)

  [s,t] = size(a);
  if s < NumChar,
    b = zeros(NumChar,NumChar);
    b(1:s,1:t) = a;
    a = b;
  end

  newText = '';
  newText = [newText sprintf('[')];
  for i = 1:NumChar,
    newText = [newText sprintf('[')];
    for j = 1:(NumChar-1),
      newText = [newText sprintf('%0.9f,', a(i,j))];
    end
    newText = [newText sprintf('%0.9f]',a(i,NumChar))];
    if i < NumChar,
      newText = [newText sprintf(',')];
    end
  end
  newText = [newText sprintf(']')];

% ---------------------------------------------------------------------------

function [newText] = subWriteGrouped(fid, a)

  newText = '';
  newText = [newText sprintf('[')];
  c = 1;
  for i = 1:4,
    newText = [newText sprintf('[')];
    for j = 1:3,
      newText = [newText sprintf('%0.9f,', a(c))];
      c = c + 1;
    end
    newText = [newText sprintf('%0.9f]',a(c))];
    c = c + 1;
    if i < 4,
      newText = [newText sprintf(',')];
    end
  end
  newText = [newText sprintf(']')];

