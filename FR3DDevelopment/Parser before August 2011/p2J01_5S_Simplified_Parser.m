% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if exist('File'),
  [File,FIndex] = zAddNTData('2j01',2,File);
else
  [File,FIndex] = zAddNTData('2j01',2);
end

F = File(FIndex);

nMin = zIndexLookup(F,'1(B)');    % start of 5S chain
nMax = zIndexLookup(F,'119(B)');  % end of 5S chain

fprintf('Original secondary structure\n');
zSecondaryStructure(F,nMin,nMax);

E = abs(fix(full(F.Edge(nMin:nMax,nMin:nMax))));
EdgeMat = (E == 1) + 2*(E > 1).*(E < 18); 
EdgeMat = triu(EdgeMat);
figure(1)
clf
caxis([0 2]);
colormap([[1 1 1]; [0 0 1]; [1 0 0]]);
pcolor(EdgeMat);
axis ij
shading flat
grid on
title('Before redefinitions');
E1 = EdgeMat;


% The 2J01 structure seems to be missing several basepairs!


fprintf('Simplified secondary structure\n');
zSecondaryStructure(F,nMin,nMax);

E = abs(fix(full(F.Edge(nMin:nMax,nMin:nMax))));
EdgeMat = (E == 1) + 2*(E > 1).*(E < 18); 
EdgeMat = triu(EdgeMat);

%fprintf('Difference %4.1f\n',max(max(EdgeMat - E1)));

figure(2)
clf
caxis([0 2]);
colormap([[1 1 1]; [0 0 1]; [1 0 0]]);
pcolor(EdgeMat);
axis ij
shading flat
grid on
title('After redefinitions');
drawnow

close all

fprintf('Nodes identified\n');
Node = pMakeNodes(F,nMin,nMax);
Node = pShiftNodeIndices(Node,nMin);

Node(1).basehalfwidth = 150;

for n=1:length(Node),
  if (n > 40) && strcmp(Node(n).type,'Hairpin'),
    Node(n).subtype = 'GNRA';
  end
end

pWriteNodesForJava(F,Node);

break

%pDisplayNodes(F,Node)

Sequence = pReadFASTA('5S_Rfam_Archaea_seed_Jesse_2_20_05.fasta');
Sequence = pSequenceConstraints(Node,Sequence);

Filename = 'Parser/5S_archaeal';

Sequence = Sequence(1:10);            % fewer sequences

disp('Parsing')
numtodo = 10;
k = 1;

for j=k:min(length(Sequence),(k+numtodo-1)),
  [Node,Sequence] = pParseSequences(Node,Sequence,j,j);
  k = j;
  Sequence = pDisplayMultipleAlignment(Node,Sequence);
  drawnow
  save(Filename,'Node','Sequence','k');
end

break

U   1
U   2
A   3
           C 122
           C 121
           A 120
G   4  cWW C 119
G   5  cWW C 118
C   6  cWW G 117
G   7  cWW C 116
G   8  cWW C 115
C   9  cWW G 114
C  10

Junction

Loop 1 - Nucleotides 11 to 68, length  58
A  11  tSs G  68
C  12
A  13
G  14  cWW C  67
C  15  cWW G  66
           A  65
G  16  cWW C  64
G  17  cWW C  63
U  18  cWW A  62
G  19  cWW C  61
G  20  cWW C  60
G  21  cWW C  59
G  22              Left  G  22  tSW C  26 | 
U  23
U  24
G  25
C  26              Left  C  26  tWS G  22 | 
C  27
           G  58
           A  57   Right A  57  tWS U  28 | 
           A  56   Right A  56  cWS C  29 | 
           U  55
U  28  cWW A  54   Left  U  28  tSW A  57 | 
C  29  cWW G  53   Left  C  29  cSW A  56 | 
           A  52   Right A  52  cWS C  30 | A  52  cHS A  51 | 
           A  51   Right A  51  cSH A  52 | 
C  30  cWW G  50   Left  C  30  cSW A  52 | 
C  31  cWW G  49
G  32  cWW C  48
U  33  cWW A  47   Right A  47  cSW C  36 | 
A  34
C  35
C  36  cWS A  47   Right A  47  cWW U  33 | 
           C  46
           A  45   Right A  45  cWS A  38 | 
           A  44   Right A  44  tHS U  39 | 
C  37  cWW G  43
A  38              Left  A  38  cSW A  45 | 
U  39              Left  U  39  tSH A  44 | 
C  40
C  41
C  42

Loop 2 - Nucleotides 69 to 114, length  46
           C 113
U  69  cWW U 112
U  70  cWW U 111
C  71  cWW G 110
C  72  cWW G 109
G  73  cWW C 108
G  74  cWW C 107
G  75  cWW C 106
G  76  tSH A 105
A  77  tHH A 104
G  78              Left  G  78  cSH U  79 | 
U  79  tWH A 103   Left  U  79  cHS G  78 | 
A  80  tHS G 102
C  81  cWW G 101
U  82  cWW G 100
G  83  cWW U  99
G  84  cWW C  98
A  85  cWW U  97
G  86  cWW C  96
U  87
G  88  cWW C  95
C  89  cWW G  94
G  90  tSH A  93
C  91
G  92

for n = 1:length(Node),
  Node(n)
  switch Node(n).type,
  case 'Basepair'
    [n Node(n).lpar(1) Node(n).rpar(1) Node(n).lip(1,1) Node(n).rip(1,1)]
  end
end

