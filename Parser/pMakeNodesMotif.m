% pMakeNodesMotif inserts a known motif into the model

if Verbose > 0,
  fprintf('Nucleotide %s%s is part of motif %s\n', File.NT(a).Base, File.NT(a).Number, File.Nucl(a).Motif.Name);
end

if Verbose > 1,
  Node(n)
  File.Nucl(a).Motif(1)
end


%File.Nucl(a).Motif(1)

ModelIndex = File.Nucl(a).Motif(1).Index;
Indices    = File.Nucl(a).Motif(1).Indices;

[y,p] = sort(Indices);               % some IL's are rotated 180 degrees

if any(diff(Indices) < 0)
  fprintf('***********************************************************\n');
  fprintf('Indices in %s are different than in %s.\n', File.Filename, File.Motifs(ModelIndex).Name);
end

NewParam = Param;
NewParam(1) = NewParam(1)-1;         % be one step less verbose

[MotifNode,Trunc] = pMakeModelFromSearchSaveFile(File.Motifs(ModelIndex).Name,NewParam);

MotifNode = MotifNode(2:end);        % omit Initial node; awkward solution
if strcmp(MotifNode(end).LeftLetter,'*'),
  MotifNode = MotifNode(1:(end-1));  % remove artificial hairpin from IL
end

for nn = 1:length(MotifNode),
  Node(n+nn) = MotifNode(nn);
  Node(n+nn).nextnode   = n+nn+1;
  Node(n+nn).LeftIndex  = y(MotifNode(nn).LeftIndex);
  Node(n+nn).RightIndex = y(MotifNode(nn).RightIndex);
  if ~isempty(Node(n+nn).MiddleIndex),
    Node(n+nn).MiddleIndex = y(Node(n+nn).MiddleIndex);
  end
  Node(n+nn).Comment = [' // Node from model: ' File.Motifs(ModelIndex).Name ' for nucleotide ' File.NT(a).Base num2str(File.NT(a).Number) ' ' Node(n+nn).Comment];

  if Verbose > 1,
    Node(n+nn)
    %pause
  end

  if Verbose > 0,
    LI = Node(n+nn).LeftIndex;
    RI = Node(n+nn).RightIndex;
    MI = Node(n+nn).MiddleIndex;

    switch Node(n+nn).type,
    case 'Initial'
      fprintf('%3d Initial   %s and %s\n', n+nn, File.NT(LI).Number, File.NT(RI).Number);

    case 'Basepair'
 
      fprintf('%3d Basepair  %4s %4s %s%s %s',n+nn, File.NT(LI).Number, File.NT(RI).Number,File.NT(LI).Base,File.NT(RI).Base,zEdgeText(File.Edge(LI,RI)));
      fprintf('In the model: %3d Basepair %s\n',nn, MotifNode(nn).Comment);

    case 'Hairpin'

      fprintf('%3d Hairpin   %4s:%4s        ',n+nn, File.NT(min(MI)).Number, File.NT(max(MI)).Number);
      fprintf('%3d Hairpin   %s\n',nn, MotifNode(nn).Comment);

      EndLoop = 1;

    end
  end
end

n = n + length(MotifNode);

if strfind(File.Motifs(ModelIndex).Name,'HL'),
  EndLoop = 1;
else
  % reset counters to move just beyond the current motif
  % because of how internal loops can be rotated 180 degrees,
  % we need to look at the ends of each strand 
  B = max(File.Nucl(a).Motif(1).Indices([(Trunc - 1) end])) - 1;
  a = min(File.Nucl(a).Motif(1).Indices([(Trunc - 1) end])) + 1;
  B = double(B);
  a = double(a);
end

