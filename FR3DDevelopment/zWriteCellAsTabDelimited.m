
% [n,t,r] = xlsread('PDB_File_Information 2009-07-24.xls');
% zWriteCellAsTabDelimited(r,'Test_2_zWriteCellAsTabDelimited.txt')

function [void] = zWriteCellAsTabDelimited(T,Filename)

fid = fopen(Filename,'w');

[a,b] = size(T);

for i = 1:a,
  for j = 1:b,

    if strcmp(class(T{i,j}),'char'),
      fprintf(fid,'%s',T{i,j});
    elseif strcmp(class(T{i,j}),'double'),
      if fix(T{i,j}) == T{i,j},
        fprintf(fid,'%d',T{i,j});
      else
        fprintf(fid,'%8.6f',T{i,j});
      end
    end


    if j < b,
      fprintf(fid,'\t');
    else
      fprintf(fid,'\n');
    end
  end
end

fclose(fid);
