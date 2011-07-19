
function [void] = zWriteMatrixAsCSV(M,FN,Decimal)



fid = fopen(FN,'w');

if Decimal == 0,
  for i = 1:length(M(:,1)),
    for j = 1:(length(M(1,:))-1),
      fprintf(fid,'%d,',M(i,j));
    end
    j = j + 1;
    fprintf(fid,'%d\n',M(i,j));
  end
else
  for i = 1:length(M(:,1)),
    for j = 1:(length(M(1,:))-1),
      fprintf(fid,'%0.4f,',M(i,j));
    end
    j = j + 1;
    fprintf(fid,'%0.4f\n',M(i,j));
  end
end

fclose(fid);


