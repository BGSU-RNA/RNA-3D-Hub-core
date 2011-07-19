% zDisplayFastaAlignment(Sequence,Columns) displays certain Columns of the
% Fasta sequences stored in Sequence.  Columns with all - entries are
% suppressed.

% zDisplayFastaAlignment(Sequence,Sequence(1).Start:Sequence(1).End)

function [void] = zDisplayFastaAlignment(Sequence,Columns)

D = zeros(1,length(Columns));

for k = 1:length(Sequence),
  D = max(D,(Sequence(k).Fasta(Columns) ~= '-'));
end

for k = 1:length(Sequence),
  fprintf('%s\n',Sequence(k).Fasta(Columns(find(D))));
end
