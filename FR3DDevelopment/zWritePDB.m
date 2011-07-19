% zWritePDB(File) writes the nucleotides in File to a PDB file

function [void] = zWritePDB(File,Filename,Rot,sh)

if nargin < 3,
  Rot = eye(3);
  sh = [0 0 0];
end

if isfield(File.NT(1),'ModelNum'),
  ModelNum = cat(2,File.NT.ModelNum);
else
  ModelNum = ones(length(File.NT),1);
end

fid = fopen(Filename,'w');       % open for writing

if max(ModelNum)-min(ModelNum) > 0,
  fprintf(fid,'MODEL   %d\n',ModelNum(1));
end

a = 1;                                         % atom number

for n = 1:length(File.NT),
  if max(ModelNum)-min(ModelNum) > 0,
    if ModelNum(n) == min(ModelNum),    % only write first model
      a = zWriteNucleotidePDB(fid,File.NT(n),a,0,Rot,sh);
    end
  else
    a = zWriteNucleotidePDB(fid,File.NT(n),a,0,Rot,sh);
  end

end

if max(ModelNum)-min(ModelNum) > 0,
  fprintf(fid,'ENDMDL\n');
end

fclose(fid);
