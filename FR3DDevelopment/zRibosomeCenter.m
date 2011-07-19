
% File = zAddNTData('EC_R1_2AVY_2AW4');

F  = zAddNTData('2avy');
FF = zAddNTData('2aw4');

% cat(2,File.NT.Chain)
% File has the 5S as chain A, 23S as chain B, 16S as chain C.

NewF.Filename = '2aw4_2avy.pdb';
for i = 1:length(F.NT),
  F.NT(i).Chain = 'C';
end
NewF.NT = [FF.NT F.NT];
zWritePDB(NewF,NewF.Filename);

C = mean(cat(1,NewF.NT.Center))

File = zAddNTData('2aw4_2avy',0,[],3)

break

VP.Sugar = 0;
VP.LabelBases = 0;
clf
zDisplayNT(File,'_A',VP);
hold on
zDisplayNT(FF,'_A',VP);

clf
zDisplayNT(File,'10:100_C',VP);
hold on
zDisplayNT(F,'10:100_A',VP);

C = mean(cat(1,File.NT.Center));
clf
[D,sh,R] = zSuperimposeNucleotides(File,'10:100_C',F,'10:100_A',VP,80);

clf
zDisplayNT(FF,[],VP);
zDisplayNT(F,[],VP);

