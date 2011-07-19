% pListSequenceFileNames looks in the Sequences folder in the current directory, finds all .fasta files, and writes their names out to a file

LoopType = 'IL';

a = dir(['Sequences' filesep '*.fasta']);

fid = fopen(['Sequences' filesep LoopType '_Sequences.txt'],'w');

for i = 1:length(a),
  b = a(i).name;
  if strcmp(a(i).name(1:2),LoopType),
    fprintf(fid,'%s\n',a(i).name);
  end
end

fclose(fid);

