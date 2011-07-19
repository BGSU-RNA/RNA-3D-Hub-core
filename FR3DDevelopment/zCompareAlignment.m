% zCompareAlignment(File,i1,i2) plots the pairwise correspondences

function [void] = zCompareAlignment(File,i1,i2,j1,j2,n1,n2)

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

View(8) = 0;                              % don't put gaps in circular diagram
View(9) = 1;                              % display the numbers

Thickness = 0.2;

clf

r = 0.8;

Matrix1 = sparse(i1,i2,ones(1,length(i1)));
Matrix1(length(File(1).NT),length(File(2).NT)) = 0;

Matrix2 = sparse(j1,j2,ones(1,length(j1)));
Matrix2(length(File(1).NT),length(File(2).NT)) = 0;

A = zNumberCircularDiagram(File(1),View,Thickness,1);
B = zNumberCircularDiagram(File(2),View,Thickness,r);
text(-1.2,1.2,File(1).Filename,'HorizontalAlignment','Left');
text(-0.2,0.2,File(2).Filename,'HorizontalAlignment','Left');

text(-1.3, 1.4, ['Comparison between ' n1 ' and ' n2]);

[g1,g2,s] = find(Matrix2 .* Matrix1 == 1);
for j = 1:length(g1),
  plot([cos(A(g1(j))) r*cos(B(g2(j)))], [sin(A(g1(j))) r*sin(B(g2(j)))],'k','LineWidth',Thickness);
  hold on
end

text(-1.3,-1.4,['Black lines indicate the ' num2str(length(g1)) ' agreements']);

[g1,g2,s] = find(Matrix2 > Matrix1);
for j = 1:length(g1),
  plot([cos(A(g1(j))) r*cos(B(g2(j)))], [sin(A(g1(j))) r*sin(B(g2(j)))],'c','LineWidth',Thickness);
  hold on
end

text(-1.3,-1.8,['Cyan lines indicate the ' num2str(length(g1)) ' extra correspondences in ' n2]);

[g1,g2,s] = find(Matrix1 > Matrix2);
for j = 1:length(g1),
  plot([cos(A(g1(j))) r*cos(B(g2(j)))], [sin(A(g1(j))) r*sin(B(g2(j)))],'r','LineWidth',Thickness);
  hold on
end

text(-1.3,-1.6,['Red lines indicate the ' num2str(length(g1)) ' missing correspondences in ' n2]);

axis equal
axis([-1.2 1.2 -2 1.2]);
axis off





