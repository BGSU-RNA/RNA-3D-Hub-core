% zHistogramDiscrepanciesInAlignment(File,i1,i2) does many local superpositions according to the alignment of File(1) and File(2) expressed by i1, i2

function [discrep] = zHistogramDiscrepanciesInAlignment(File,i1,i2)

itoa     = zeros(1,length(File(1).NT));
itoa(i1) = 1:length(i1);                 % maps indices in File(1) to alignment

discrep = [];

for j = 1:length(i1),
  d = File(1).Distance(i1(j),:);
  k = find((d > 0) .* (d < 8));
  m = itoa([i1(j) k]);                   % entries of i1
  m = m(find(m));                        % only keep non-zero entries
  if length(m) > 1,
    discrep(j) = xDiscrepancy(File(1),i1(m),File(2),i2(m));
  end
end

