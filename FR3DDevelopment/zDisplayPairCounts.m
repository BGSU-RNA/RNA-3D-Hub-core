% zDisplayPairCounts displays counts of pairwise interactions, either from exemplars or from structures.  The format is Count(base1code,base2code,category).

function [T] = zDisplayPairCounts(Count,CountFilename,Letters)

if nargin < 3,
  Letters = 'ACGU';
end

L = length(Letters);


Total = [];

[s,t,MaxCat] = size(Count);

if MaxCat > 100,
  m = 10;                           % category numbers are all multiplied by 10
else
  m = 1;
end

for ca = 1:12,
  Total(ca) = sum(sum(Count(:,:,m*ca)));  % total for each basepair family
end

row = 1;
clear T

% ---------------------------------------- counts by family

T{row,1} = 'Family';
T{row,2} = 'Count';
T{row,3} = '% of total';

for ca = 1:12,                           % families 1-12 only
  T{row+ca,1} = zEdgeText(ca,0);
  T{row+ca,2} = Total(ca);
  T{row+ca,3} = 100*Total(ca)/sum(Total(1:12));
end

row = row + ca + 2;

S = 6;

% ---------------------------------------- counts by basepair within family
%                                          L lines for each family

for ca = 1:12,
  T{row,1} = zEdgeText(ca,0);
  T{row,1+S} = zEdgeText(ca,0);
  for i = 1:L,
    T{row+i,1} = Letters(i);
    T{row,1+i} = Letters(i);
    for j = 1:L,
      T{row+i,1+j} = Count(i,j,m*ca);
    end

    T{row+i,1+S} = Letters(i);
    T{row,1+i+S} = Letters(i);
    for j = 1:L,
      T{row+i,1+j+S} = 100*Count(i,j,m*ca)/Total(ca);
    end

  end

  row = row + 6;
end

% ---------------------------------------- counts by basepair within family
%                                          one line for each family

T{row,2} = 'Counts from 3D structures';
row = row + 1;
T{row,18} = 'Total';
T{row,19} = '% of total';
for i = 1:L,
  for j = 1:L,
    pc = L*(j-1) + i;
    T{row,1+pc} = [Letters(j) Letters(i)];
  end
end

for ca = 1:12,
  T{row+ca,1} = zEdgeText(ca,0);
  for i = 1:L,
    for j = 1:L,
      pc = L*(j-1) + i;
      if Count(i,j,m*ca) ~= 0,
        T{row+ca,1+pc} = Count(j,i,m*ca);
      end
    end
  end
  T{row+ca,18} = Total(ca);
  T{row+ca,19} = 100*Total(ca)/sum(Total(1:12));
end

row = row + 15;

% --------------------------------------- frequencies by basepair within family
%                                         one line for each family

T{row,2} = 'Frequencies from 3D structures';
row = row + 1;
T{row,18} = 'Total';
T{row,19} = '% of total';
for i = 1:L,
  for j = 1:L,
    pc = L*(j-1) + i;
    T{row,1+pc} = [Letters(j) Letters(i)];
  end
end

for ca = 1:12,
  T{row+ca,1} = zEdgeText(ca,0);
  for i = 1:L,
    for j = 1:L,
      pc = L*(j-1) + i;
      if Count(i,j,m*ca) ~= 0,
        T{row+ca,1+pc} = 100*Count(j,i,m*ca)/Total(ca);
      end
    end
  end
  T{row+ca,18} = 100;
  T{row+ca,19} = 100*Total(ca)/sum(Total(1:12));
end

row = row + 15;

% --------------------------------------- counts by basepair

T{row,2} = 'Counts by pair, from 3D structures';
row = row + 1;

for ca = 1:12,
  T{row,ca+1} = zEdgeText(ca,0);
end
T{row,14} = 'Total';
T{row,15} = '% of total';

for i = 1:L,
  for j = 1:L,
    pc = L*(i-1) + j;
    T{row+pc,1} = [Letters(i) Letters(j)];
    T{row+pc,14} = sum(Count(i,j,m*(1:12)));
    T{row+pc,15} = 100*sum(Count(i,j,m*(1:12)))/sum(Total(1:12));
    for ca = 1:12,
      T{row+pc,ca+1} = Count(i,j,m*ca);
    end
  end
end

row = row + 19;

% --------------------------------------- counts by basepair

T{row,2} = 'Frequencies with which each pair falls in each family, from 3D structures';
row = row + 1;

for ca = 1:12,
  T{row,ca+1} = zEdgeText(ca,0);
end
T{row,14} = 'Total';

for i = 1:L,
  for j = 1:L,
    pc = L*(i-1) + j;
    T{row+pc,1} = [Letters(i) Letters(j)];
    total = sum(Count(i,j,m*(1:12)));
    for ca = 1:12,
      T{row+pc,ca+1} = 100*Count(i,j,m*ca)/total;
    end
    T{row+pc,14} = 100;
  end
end

row = row + 19;

% --------------------------------------- counts by basepair

T{row,2} = 'Counts by pair, from 3D structures';
row = row + 1;

for ca = 1:12,
  T{row,ca+1} = zEdgeText(ca,0);
end
T{row,14} = 'Total';
T{row,15} = '% of total';

k = 0;

for i = 1:L,
  for j = i:L,
    k = k + 1;
    s = 0;
    if i < j,                                  % in alphabetical order
      for ca = 1:12,
        T{row+k,ca+1} = Count(i,j,m*ca) + Count(j,i,m*ca);
        s = s + T{row+k,ca+1};
      end
      T{row+k,1} = [Letters([i j]) '/' Letters([j i])];
    else
      for ca = 1:12,
        T{row+k,ca+1} = Count(i,j,m*ca);
        s = s + T{row+k,ca+1};
      end
      T{row+k,1} = Letters([i j]);
    end
    T{row+k,14} = s;
    T{row+k,15} = 100*s/sum(Total(1:12));
  end
end

row = row + 13;

% --------------------------------------- counts by basepair

T{row,2} = 'Frequency with which each pair falls in each family, from 3D structures';
row = row + 1;

for ca = 1:12,
  T{row,ca+1} = zEdgeText(ca,0);
end
T{row,14} = 'Total';

k = 0;

for i = 1:L,
  for j = i:L,
    k = k + 1;
    s = 0;
    if i < j,                                  % in alphabetical order
      total = sum(Count(i,j,m*(1:12))) + sum(Count(j,i,m*(1:12)));
      for ca = 1:12,
        T{row+k,ca+1} = 100*(Count(i,j,m*ca) + Count(j,i,m*ca))/total;
      end
      T{row+k,1} = [Letters([i j]) '/' Letters([j i])];
    else
      total = sum(Count(i,j,m*(1:12)));
      for ca = 1:12,
        T{row+k,ca+1} = 100*Count(i,j,m*ca)/total;
      end
      T{row+k,1} = Letters([i j]);
    end
    T{row+k,14} = 100;
  end
end

% --------------------------------------- write to spreadsheet

if ~isempty(CountFilename),
  xlswrite(CountFilename,T);
end


