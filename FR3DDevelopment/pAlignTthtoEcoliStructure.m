
File1 = zAddNTData({'1j5e'});
File2 = zAddNTData({'2avy'});

S1 = cat(2,File(1).NT.Base);                % extract sequence from file
S2 = cat(2,File(2).NT.Base);                % extract sequence from file

% in the future, extract the sequences from here and run the Java program

% for now, use the sequences pasted into this program:

pInferVsActual3DStructure

% infer the indices aligned from the sequence alignment:

[i1,i2] = pInferAlignmentFromSequenceAlignment(File1,ES,File2,JS);

% write to an Excel spreadsheet:

rWriteAlignmentBPComparison(File1,1:length(File1.NT),File2,1:length(File2.NT),i1,i2)

