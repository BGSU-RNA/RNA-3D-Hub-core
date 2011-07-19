% pParse 
warning off
Filename = '5S_June_Euryarchaeota_B.mat';

k = 1;
if exist(Filename,'file') > 0,            % load previous work, if it exists
  load(Filename);
  k = k + 1;                              % go to next sequence
  fprintf('Loaded previous computations\n');
else                                      % define model, load sequences
  clear Node
%  Sequence = pReadFASTA('5S_Rfam_Full.fasta');
  Sequence = pReadFASTA('RFAMfull_Euryarchaeota.fasta');
%  Sequence = Sequence(1:1000);
  Nodes_5SA_June
  [Node,Sequence]=pMakeConstraints(Node,Sequence);
  disp('Nodes Created')
end

for j=1:length(Sequence),
%  Sequence(j).Fasta = '';
end

numtodo = 5000;

for j=k:min(length(Sequence),(k+numtodo-1)),
  [Node,Sequence] = pParseSequences(Node,Sequence,j,j,0);
  k = j;
  drawnow

  if mod(k,10) == 0,
    disp('Saving Parse')
    disp(k)
    save(Filename,'Node','Sequence','k');
    pDisplayMultipleAlignment(Node,Sequence(max(1,j-30):j));
  end
end

%Sequence=pDisplayMultipleAlignment(Node,Sequence);

return

OK = [];
for j=1:k,
  if isempty(Sequence(j).TraceInfo),
    OK(j) = 0;
  else
    OK(j) = 1;
  end 
end

m = find(OK);

pDisplayMultipleAlignment(Node,Sequence(m));
pWriteFasta(Sequence(m), 'RFAM_alignment_of_parsed_archaeal_sequences.fasta');
