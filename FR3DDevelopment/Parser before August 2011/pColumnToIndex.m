% pColumnToIndex maps column of an alignment to index of each letter

function [cti,itc] = pColumnToIndex(S)

cti = zeros(1,length(S));                % column to index map
i = 1;

for c = 1:length(S),
  cti(c) = i;
  if ismember(S(c),'ACGUacgu'),
    itc(i) = c;
    i = i + 1;
  end
end
