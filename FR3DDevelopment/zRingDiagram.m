% zRingDiagram(File,i1,i2,c) plots the alignment of File(1) and File(2) indicated by vectors i1 and i2, coloring each line segment according to c

function [void] = zRingDiagram(File,i1,i2,c)

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

View(8) = 0;                              % don't put gaps in circular diagram
View(9) = 1;                              % display the numbers

Thickness = 0.2;                          

clf

r = 0.8;                                  % inner radius

A = zNumberCircularDiagram(File(1),View,Thickness,1);
B = zNumberCircularDiagram(File(2),View,Thickness,r);
text(-1.2,1.2,File(1).Filename,'HorizontalAlignment','Left');
text(-0.2,0.2,File(2).Filename,'HorizontalAlignment','Left');

% text(-1.3, 1.4, ['Comparison between ' n1 ' and ' n2]);

[s,t] = size(c);

colo = c;

if s == 1,
  for j = 1:length(i1),
    colo(j,1) = c;
  end
end

if (t == 1) && strcmp(class(c),'double'),
  colormap('default');
  map = colormap;
  mins = min(c);
  maxs = max(c);
  colo = map(7+ceil(48*(c-mins)/(maxs-mins)),:);
end

  

for j = 1:length(i1),
  plot([cos(A(i1(j))) r*cos(B(i2(j)))], [sin(A(i1(j))) r*sin(B(i2(j)))],'Color',colo(j,:),'LineWidth',Thickness);
  hold on
end

axis equal
axis([-1.2 1.2 -2 1.2]);
axis off





