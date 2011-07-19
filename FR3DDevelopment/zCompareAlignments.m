% zCompareAlignments(Alignment1,Alignment2,File) does simple comparisons
% between alignments Alignment1 and Alignment2 of the two files in File

% An Alignment is an M element cell array of vectors of nucleotide indices.
% Each element indicates one, two, or more indices of each file that align

function [void] = zCompareAlignments(Alignment1,Alignment2,File)

AMatrix1 = sparse(zeros(File(1).NumNT,File(2).NumNT));

for i = 1:length(Alignment1.Correspondence),
  for j = 1:length(Alignment1.Correspondence(i).File1),
%[i j]
    AMatrix1(Alignment1.Correspondence(i).File1(j),Alignment1.Correspondence(i).File2(j)) = 1;
  end
end

AMatrix2 = sparse(zeros(File(1).NumNT,File(2).NumNT));

for i = 1:length(Alignment2.Correspondence),
  for j = 1:length(Alignment2.Correspondence(i).File1),
    AMatrix2(Alignment2.Correspondence(i).File1(j),Alignment2.Correspondence(i).File2(j)) = 1;
  end
end

% --------------------- Diagnostic information comparing the two alignments

figure(1)
clf
spy( AMatrix1 .* (AMatrix1 == AMatrix2));
title('Locations of agreement between the two alignments');

figure(2)
clf
spy( AMatrix1 .* (AMatrix1 ~= AMatrix2));
title('Elements of Alignment1 that are not in Alignment2');

figure(3)
clf
spy( AMatrix2 .* (AMatrix1 ~= AMatrix2));
title('Elements of Alignment2 that are not in Alignment1');

