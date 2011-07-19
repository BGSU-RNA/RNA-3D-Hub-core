% zMakeScoreMatrix uses 3D superposition to make a score matrix for
% aligning nucleotides from two files

% File = zAddNTData({'1s72', '2j01'},2);
% Indices1 = 2750:2871;
% zMakeScoreMatrix(File,Indices1,2773:2891);

function [S] = zMakeScoreMatrix(File,Indices1,Indices2)

m = 6;                  % number of nucleotides to superimpose at once

N1 = length(Indices1);
N2 = length(Indices2);

% find indices of nucleotides near each nucleotide

for k = 1:N1,
  a     = find(File(1).Distance(Indices1(k),:)); % nucleotides close to this NT
  [y,i] = sort(File(1).Distance(Indices1(k),a));
  b     = [Indices1(k) a(i(1:min(m-1,length(i))))];
  I1{k} = sort(b);
-
%clf
%zDisplayNT(File(1),I1{k});
%pause

end

for k = 1:N2,
  a     = find(File(2).Distance(Indices2(k),:)); % nucleotides close to this NT
  [y,i] = sort(File(2).Distance(Indices2(k),a));
  b     = [Indices2(k) a(i(1:min(m-1,length(i))))];
  I2{k} = sort(b);
end

VP.Plot  = 0;
VP.Write = 0;

for k = 1:N1,
  for j = 1:N2,

%if abs(k-j) < 3,
    S(k,j) = zSuperimposeNucleotides(File(1),I1{k},File(2),I2{j},VP);

%fprintf('Discrepancy %6.3f\n', S(k,j));
%drawnow
%pause
%end

  end
end




