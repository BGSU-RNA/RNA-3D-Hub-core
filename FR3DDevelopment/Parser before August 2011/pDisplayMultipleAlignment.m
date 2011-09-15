% pDisplayMultipleAlignment displays the multiple alignment of parsed sequences

function [Sequence] = pDisplayMultipleAlignment(Node,Sequence,Group)

N = length(Node);                         % number of nodes
S = length(Sequence);

if nargin == 2,
  Group = ones(1,S);
end

% ------------------------------------- Determine space to leave for insertions

for n=1:N,
  maxinsert(n).left    = 0;               % largest insertion on left, node n
  maxinsert(n).right   = 0;               % largest insertion on right, node n
  maxinsert(n).middle  = 0;
  maxinsert(n).hairpin = 1;               % longest hairpin
end



for k=1:S,                                 % loop through all parsed sequences
 if ~isempty(Sequence(k).TraceInfo),
  x = Sequence(k).X;                       % kth sequence
  L = length(x);
%k
  maxinsert = pFindMaxInsertions(Node,Sequence(k).TraceInfo,maxinsert,L);
 end
end
save('MaxInsert.mat','maxinsert');
% for k=30:30
%   k
%   aL=maxinsert(k).left;
%   aR=maxinsert(k).right;
%   aM=maxinsert(k).middle;
%   aH=maxinsert(k).hairpin;
%   aL
%   aR
%   aM
%   aH
%   disp('--------------')
% end

% ------------------------------------- Display header for multiple alignment

header = pMakeLoopHeader(Node,maxinsert,1);
fprintf('     %s\n',header{3});
fprintf('     %s\n',header{2});
fprintf('     %3s  Log prob  Organism\n',header{1});

% ------------------------------------- Display sequences aligned
Sequence(1).header=header;
for k=1:S,                                % loop through sequences
drawnow
 if ~isempty(Sequence(k).TraceInfo),
  x = Sequence(k).X;
  L = length(x);

  n = 1;
  s = 1;                                  % state to start in
  i = 1;                                  % left entry of x we're explaining
  j = L;                                  % right entry of x we're explaining
%   k
  alignment = pAlignSequence(Node,Sequence(k).TraceInfo,maxinsert,n,s,i,j,x);
  Sequence(k).Alignment = alignment;

  fprintf('%4s %s %s %8.2f  %s\n',Sequence(k).FastaNum,alignment,Sequence(k).FastaNum,Sequence(k).TraceInfo(1,1).mp,Sequence(k).Organism);
   if 1-strcmp(alignment(find(alignment~='-' & alignment~='*' & alignment~='+')),Sequence(k).X)
     fprintf('Error displaying sequence %s (index %d)\n',Sequence(k).FastaNum,k);
     fprintf('%s\n',alignment(find(alignment~='-' & alignment~='*' & alignment~='+')));
     XX=Sequence(k).X;
     fprintf('%s\n',XX);
  end 
  if Group(k) ~= Group(min(S,k+1)),
    fprintf('\n');
  end
  

 end
end
