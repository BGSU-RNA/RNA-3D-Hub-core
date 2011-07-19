% zCircularDiagram(File,Thickness,View) plots the pairwise interactions in File
% using colored chords around a unit circle.  Thickness controls the
% thickness of the lines, for different graphic output formats.

% zCircularDiagram('1s72') will load the file and display the diagram.

% Helpful suggestions for how to save the figure as a graphic file:
%  clf
%  zCircularDiagram(File(f),1);
%  saveas(gcf,[mypath FN '_circular_diagram.png'],'png');
%  [X,map] = imread([mypath FN '_circular_diagram.png']);
%  Y = X(30:830,210:1030,:);
%  imwrite(Y,[mypath FN '_circular_diagram.png']);

%  clf
%  zCircularDiagram(File(f),0.1);
%  saveas(gcf,[mypath FN '_circular_diagram.pdf'],'pdf');

% clf
% zCircularDiagram('1s72_9',1,[1 1 1 1 1 1 1 0 1 0 0])
% saveas(gcf,'1s72_9_circular_diagram.png','png');

%  zCircularDiagram(File,1,[1 1 0 0 0 1 1]);

%  View tells what to display
%  1 = nested cWW
%  2 = nested non-cWW
%  3 = non-nested cWW
%  4 = non-nested non-cWW
%  5 = stack
%  6 = BPh
%  7 = explanatory text
%  8 = leave gaps in diagram for breaks in nucleotide chain
%  9 = display nucleotide numbers in the diagram
% 10 = display nucleotide letters in the diagram
% 11 = display B for BPh interactions in the diagram

function [Tally] = zCircularDiagram(File,Thickness,View)

Tally = zeros(1,6);

if nargin < 2,
  Thickness = 1;
end

if nargin < 3,
  View = [1 1 1 1 1 1 1 1 1 1 1];
end

while length(View) < 9,
  View = [View 1];
end

if strcmp(class(File),'char'),
  a = strfind(File,'_');
  if ~isempty(a),
    OrigFilename = strrep(File,'_','\_');
    chain = File([a(1) a(1)+1]);
    Filename = File(1:(a(1)-1));
    File = zGetNTData(Filename,0);
    i    = zIndexLookup(File,chain);
    File = zSubFile(File,i);
  else
    Filename = File;
    OrigFilename = File;
    File = zGetNTData(Filename,0);
  end
else
  OrigFilename = File.Filename;
end

E  = fix(abs(File.Edge));
C  = File.Crossing;

Color = zKeepZeros(zSparseValues(E,1,1),C);      % nested cWW
Color = Color + 2*zKeepZeros(zSparseRange(E,2,13,1),C); % nested non-cWW
Color = Color + 3*(E==1).*(C>0);                 % non-nested cWW
Color = Color + 4*zSparseRange(E,2,13,1).*(C>0); % non-nested non-cWW
Color = Color + 5*zSparseRange(E,21,24,1);       % stacks

clear E C

[i,j,c] = find(triu(Color));            % keep one of each interaction

BP = abs(File.BasePhosphate);           % 
BP = zSparseRange(BP,0,99,1);           % extract BPh interactions
[ii,jj,cc] = find(6*BP);                % base-phosphate interactions

k = find(ii ~= jj);                     % eliminate self BPh interactions

i = [i; ii(k)];                         % append base-phosphate interactions
j = [j; jj(k)];
c = [c; cc(k)];

BR = abs(File.BaseRibose);           % 
BR = zSparseRange(BR,0,99,1);           % extract BR interactions
[ii,jj,cc] = find(7*BR);                % base-ribose interactions

k = find(ii ~= jj);                     % eliminate self BR interactions

i = [i; ii(k)];                         % append base-ribose interactions
j = [j; jj(k)];
c = [c; cc(k)];

if length(i) > -10,
  [A,mA] = zNumberCircularDiagram(File,View,Thickness,1);

% ---------------------------------------- Draw the interactions in color

[y,k] = sort(-abs(c));               % sort by decreasing color, for overlap

i = i(k);
j = j(k);
c = c(k);

color = 'bcrgymo';

shift = 2*pi*[0 0 0 0 -0.15 0.15 0.25]/max(2000,length(File.NT));    % shift so arcs don't overlap

for k = 1:length(i),
  if View(c(k)) > 0,
    sh = shift(c(k));
    zUnitCircleArc([cos(A(i(k))+sh) cos(A(j(k))+sh)], [sin(A(i(k))+sh) sin(A(j(k))+sh)],c(k),Thickness);

%  if mod(i(k),5) == 2,
      thetai = A(i(k));
      thetaj = A(j(k));

% make the radius depend on how many have already been written here; no overlap

    r = 1.005 + 0.01*c(k);
%    text(r*cos(thetai), r*sin(thetai), num2str(File.Crossing(i(k),j(k))),'FontSize',1, 'Rotation', 0, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color',color(c(k)));
%    text(r*cos(thetaj), r*sin(thetaj), num2str(File.Crossing(i(k),j(k))),'FontSize',1, 'Rotation', 0, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color',color(c(k)));

%[i(k) j(k) thetai thetaj]

%  end

    if c(k) == 6,                 % label base-phosphate interactions
      thetai = A(i(k));
      thetaj = A(j(k));

      if cos(thetai) > 0,
        angle = 180*thetai/pi;
        ha    = 'left';
      else
        angle = 180*(thetai - pi)/pi;
        ha    = 'right';
      end

      if View(11) > 0 && Thickness < 1,
        if File.BasePhosphate(i(k),j(k)) > 0,
          text(1.03*cos(thetai), 1.03*sin(thetai), 'B','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color','m');
%          text(1.03*cos(thetaj), 1.03*sin(thetaj), 'Ph','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
        else
          text(1.03*cos(thetai), 1.03*sin(thetai), 'Ph','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
%         text(1.03*cos(thetai), 1.03*sin(thetai), 'B','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
        end
      end
    end
  end
end
axis equal
axis([-1.2 1.2 -2 1.2]);
axis off


% -------------------------- Write data out to a file

if 0 > 1,
  fid = fopen([File.Filename '_nucleotides_angles.txt'],'w');
  for w = 1:length(File.NT),
    fprintf(fid,'%d\t%7.6f\t%s\n', w, A(w), File.NT(w).Number);
  end
  fclose(fid);

  fid = fopen([File.Filename '_interactions.txt'],'w');
  for w = 1:length(i),
    fprintf(fid,'%d\t%d\t%d\n', i(w), j(w), c(w));
  end
  fclose(fid);

end




text(-1.2,1.2,OrigFilename,'HorizontalAlignment','Left');

  cww = length(find(c == 1));
  Tally(1,1) = cww;
  noncww = length(find(c == 2));
  Tally(1,2) = noncww;
  nonnestcww = length(find(c == 3));
  Tally(1,3) = nonnestcww;
  nonnestnoncww = length(find(c == 4));
  Tally(1,4) = nonnestnoncww;
  stack = length(find(c == 5));
  Tally(1,5) = stack;
  bph = length(find(c == 6));
  Tally(1,6) = bph;
  br = length(find(c == 7));
  Tally(1,7) = br;

if View(7) > 0,

  if View(1) > 0,
    text(-1.3,-1.4,['Dark blue chords indicate the ' num2str(cww) ' nested Watson-Crick basepairs']);
  end

  if View(2) > 0,
    text(-1.3,-1.6,['Cyan chords indicate the ' num2str(noncww) ' nested non-Watson-Crick basepairs']);
  end

  if View(3) > 0,
    text(-1.3,-1.5,['Red chords indicate the ' num2str(nonnestcww) ' non-nested Watson-Crick basepairs']);
  end

  if View(4) > 0,
    text(-1.3,-1.7,['Green chords indicate the ' num2str(nonnestnoncww) ' non-nested non-Watson-Crick basepairs']);
  end

  if View(5) > 0,
    text(-1.3,-1.8,['Yellow chords indicate the ' num2str(stack) ' stacking interactions']);
  end

  if View(6) > 0,
    text(-1.3,-1.9,['Magenta chords indicate the ' num2str(bph) ' base-phosphate interactions']);
  end

  if View(7) > 0,
    text(-1.3,-2.0,['Orange chords indicate the ' num2str(br) ' base-ribose interactions']);
  end
end

end
