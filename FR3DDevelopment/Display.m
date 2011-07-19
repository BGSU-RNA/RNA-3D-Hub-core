
function [void] = Display(F1,F2,F3,F4)

VP.LabelBases = 0;

close all
figure(1)
clf
F1 = zAddNTData(F1);
fprintf('%s\n',cat(2,F1.NT.Base));
zDisplayNT(F1,[],VP)


figure(2)
clf
F2 = zAddNTData(F2);
fprintf('%s\n',cat(2,F2.NT.Base));
zDisplayNT(F2,[],VP)

if nargin > 2,
  figure(3)
  clf
  zDisplayNT(F3,[],VP)
end

if nargin > 3,
  figure(4)
  clf
  zDisplayNT(F4,[],VP)
end
