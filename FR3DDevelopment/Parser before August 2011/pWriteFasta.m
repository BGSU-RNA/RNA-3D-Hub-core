
function [void] = pWriteFasta(Sequence,Filename)

fidOUT = fopen(Filename,'w+');

for k=1:length(Sequence),
  fprintf(fidOUT,'> %s\n',Sequence(k).Organism);
  c = 1;
  L = length(Sequence(k).Fasta);
  while c <= L,
    fprintf(fidOUT,'%s\n', Sequence(k).Fasta(c:min(L,c+59)));
    c = c + 60;
  end
end

fclose(fidOUT);