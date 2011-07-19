% pGetAlignedIndices(S1,S2) determines the indices of sequences S1 and S2 which are aligned to one another as inferred by characters one above the other

% ------------------------------------------- Determine alignment from JAR3D

function [i1,i2] = pGetAlignedIndices(S1,S2)

S1 = upper(S1);
S2 = upper(S2);

i1 = [];
i2 = [];

if length(S1) ~= length(S2),
  fprintf('The sequences have different lengths; cannot trust that they imply an alignment\n');
else
  a = 1;
  b = 1;
  for c = 1:length(S1),
    if ismember(S1(c),'ACGU') && ismember(S2(c),'ACGU'),
      i1 = [i1 a];
      i2 = [i2 b];
      a = a + 1;
      b = b + 1;
    elseif ismember(S1(c),'ACGUacgu'),
      a = a + 1;
    elseif ismember(S2(c),'ACGUacgu'),
      b = b + 1;
    end
  end
end
