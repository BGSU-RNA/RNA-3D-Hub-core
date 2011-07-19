% pInferAlignmentFromSequenceAlignment uses an alignment of sequences S1 and S2, which are the nucleotide sequences from File1 and File2, respectively, to infer an alignment of the nucleotides from File1 and File2.  The output is a list of indices of aligned nucleotides from the two structures.

% the program assumes that every index in File1 and File2 has an A, C, G, or U base.

function [i1,i2] = pInferAlignmentFromSequenceAlignment(S1,S2)

[cti1,itc1] = pColumnToIndex(S1);
[cti2,itc2] = pColumnToIndex(S2);

i1 = [];
i2 = [];

for c = 1:length(S1),                      % go through columns of alignment
  if any(S1(c) == 'ACGU') && (any(S2(c) == 'ACGU')),   % NTs aligned
    i1 = [i1 cti1(c)];
    i2 = [i2 cti2(c)];
  end
end
