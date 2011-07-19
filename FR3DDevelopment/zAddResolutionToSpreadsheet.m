% zAddResolutionToSpreadsheet reads a spreadsheet, looks up resolution and number of basepairs per nucleotide, and writes them back out to a new spreadsheet

[nn,tt,rr] = xlsread('temp.xlsx');

firstline = 3;
PDBIDCol  = 4;

load PDBInfo

clear res

for i = firstline:length(rr(:,1)),
  FN = rr{i,PDBIDCol};
  if ~isnan(FN),
    j = find(ismember(t(:,1),FN));
    if ~isempty(j),
      rr{i,PDBIDCol+1} = n(j(1),1);             % resolution
      rr{i,PDBIDCol+2} = n(j(1),3)/n(j(1),2);      % basepairs per nucleotide

      res(i) = n(j(1),1);
    end
  end
end

xlswrite('temp4.xlsx',rr)

hist(res(res>0),20)
title('Resolution over coplanar triples with fewer than 2 basepairs')
ax = axis;
ax(1) = 0;
axis(ax);
%saveas(gcf,'Resolution over coplanar triples with fewer than 2 basepairs.png');
mean(res(res>0))

break

allres = res;

for i = 1:length(res),
  j = find(allres == res(i));
  j = j(1);
  allres = allres([1:(j-1) (j+1):end]);
end

res = allres;

hist(res(res>0),20)
title('Resolution over coplanar triples with 2 or more basepairs')
ax = axis;
ax(1) = 0;
axis(ax);
saveas(gcf,'Resolution over coplanar triples with 2 or more basepairs.png');
mean(res(res>0))
