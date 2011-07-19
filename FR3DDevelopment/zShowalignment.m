% zShowAlignment(File,i1,i2) plots the pairwise correspondences

function [void] = zShowAlignment(File,i1,i2)

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

View(8) = 0;
Thickness = 0.5;

clf

r = 0.8;

A = zNumberCircularDiagram(File(1),View,Thickness,1);
B = zNumberCircularDiagram(File(2),View,Thickness,r);

for j = 1:length(i1),
  plot([cos(A(i1(j))) r*cos(B(i2(j)))], [sin(A(i1(j))) r*sin(B(i2(j)))],'k');
  hold on
end

text(-1.2,1.2,File(1).Filename,'HorizontalAlignment','Left');
text(-0.2,0.2,File(2).Filename,'HorizontalAlignment','Left');

