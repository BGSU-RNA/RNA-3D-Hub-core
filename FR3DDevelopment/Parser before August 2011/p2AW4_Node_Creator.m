% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if exist('File'),
  [File,FIndex] = zAddNTData('2aw4',0,File);
else
  [File,FIndex] = zAddNTData('2aw4',0);
end

F = File(FIndex);

%zSecondaryStructure(F,1,117);

E = abs(fix(full(F.Edge(1:117,1:117))));
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

F = pModifyEdge(F, '36','A', '49','A','cWS');
F = pModifyEdge(F, '56','A', '59','A','tSW');
F = pModifyEdge(F, '29','A', '58','A','cSW');
F = pModifyEdge(F, '30','A', '57','A','cSW');
F = pModifyEdge(F, '75','A','101','A','cWW');
F = pModifyEdge(F,'102','A', '74','A','bif');
F = pModifyEdge(F, '76','A','100','A','bif');
%F = pModifyEdge(F, '13','A', '69','A','tSS');
F = pModifyEdge(F, '69','A','107','A','tHS');
%F = pModifyEdge(F,'108','A', '67','A','tWH');
F = pModifyEdge(F, '66','A', '68','A','cSH');
F = pModifyEdge(F, '53','A', '52','A','cHS');

zSecondaryStructure(F,1,117);

E = abs(fix(full(F.Edge(1:117,1:117))));
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

Node = pMakeNodes(F,1,117);
%pDisplayNodes(F,Node)

Sequence = pReadFASTA('5S_Rfam_Bacteria_Zirbel.fasta');
Sequence = pSequenceConstraints(Node,Sequence);

Filename = 'Parser/5S_bacterial';

Sequence = Sequence(1:10);            % fewer sequences

k = 1;

disp('Parsing')
numtodo = 10;

for j=k:min(length(Sequence),(k+numtodo-1)),
  [Node,Sequence] = pParseSequences(Node,Sequence,j,j);
  k = j;
  Sequence = pDisplayMultipleAlignment(Node,Sequence);
  drawnow
  save(Filename,'Node','Sequence','k');
end

break

G   2  cWW C 118
C   3  cWW G 117
C   4  cWW G 116
U   5  cWW A 115
G   6  cWW C 114
G   7  cWW C 113
C   8  cWW G 112
G   9  cWW U 111
G  10  cWW C 110
C  11  cSW A 109
C  12

Junction

Loop 1 - Nucleotides 13 to 69, length  57
G  13  tSS G  69   Right G  69  tHS G 107 | 
U  14
A  15
G  16  cWW C  68   Right C  68  cHS A  66 | 
C  17  cWW G  67   Right G  67  tHW A 108 | 
           A  66   Right A  66  cSH C  68 | 
G  18  cWW U  65
C  19  cWW G  64
G  20  cWW C  63
G  21  cWW C  62
U  22  cWW G  61
G  23  cWW C  60
G  24              Left  G  24  tSW C  27 | 
U  25
C  26
C  27              Left  C  27  tWS G  24 | 
           A  59   Right A  59  tWS G  56 | 
           A  58   Right A  58  cWS A  29 | 
           A  57   Right A  57  cWS C  30 | 
C  28  cWW G  56   Right G  56  tSW A  59 | 
A  29  cWW U  55   Left  A  29  cSW A  58 | 
C  30  cWW G  54   Left  C  30  cSW A  57 | 
           A  53   Right A  53  cWS C  31 | A  53  cHS A  52 | 
           A  52   Right A  52  cSH A  53 | 
C  31  cWW G  51   Left  C  31  cSW A  53 | 
U  32  cWW A  50
G  33  cWW C  49   Right C  49  cSW C  36 | 
A  34  cWW U  48
C  35
C  36              Left  C  36  cWS C  49 | 
C  37
           C  47   Right C  47  tWS G  44 | 
           A  46   Right A  46  cWS A  39 | 
           A  45   Right A  45  tHS U  40 | 
C  38  cWW G  44   Right G  44  tSW C  47 | 
A  39              Left  A  39  cSW A  46 | 
U  40              Left  U  40  tSH A  45 | 
G  41
C  42
C  43

Loop 2 - Nucleotides 70 to 110, length  41
           A 109   Right A 109  cWS C  11 | 
           A 108   Right A 108  tWH G  67 | 
           G 107   Right G 107  tSH G  69 | 
C  70  cWW G 106
C  71  cWW G 105
G  72  tSH A 104
A  73  tHW U 103
U  74  bif G 102
G  75  cWW A 101
G  76  bif G 100
U  77  tWH A  99
A  78  tHS G  98
G  79  cWW C  97
U  80  cWW G  96
G  81  cWW U  95
U  82  cWW A  94
G  83  cWW C  93
G  84  cWW C  92
G  85  cWW C  91
G  86  cWW C  90
U  87
C  88
U  89
