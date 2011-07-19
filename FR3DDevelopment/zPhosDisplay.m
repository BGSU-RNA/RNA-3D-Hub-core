% zPhosDisplay displays base-phosphate interaction parameters in a variety of ways

function [void] = zPhosDisplay(D,Limits)

if nargin < 2,
  Limits = 0;
end

zStandardBases

L = {'A','C','G','U'};

% ------------------------------ specify cutoffs for classification

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 4.5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 4.0;                           % near category

AL = 130;                                      % angle limit for BP
nAL = 110;                                     % angle limit for nBP

DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;

nDL([4 6 11]) = nNitrogenDist;
nDL([7 8 9])  = nCarbonDist;

% ------- display interacting oxygens and phosphorus atoms together with bases

figure(1);
clf

for v = 1:4,
%  figure(v)
%  clf
  subplot(2,2,v);

%  r = find(D(:,4) == v);            % select the current nucleotide code

  r = find((D(:,4) == v) .* (abs(D(:,12)) < 3)); % select nucl and vert displ

  DD = D(r,:);                       % use only these lines of data

  r = find(DD(:,17) == 1);           % only display the best oxygen for each H
  DD = DD(r,:);

  s = 6*ones(size(DD(:,1)));        % all have the same size
%  s = 10*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
%  s = 10*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
%  s = 10*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

  scatter3(DD(:,13),DD(:,14),DD(:,15), 0.2*s, 'k', 'filled'); % phosphorus
  hold on

if Limits == 0,
  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   DD(:,9), 'filled'); % color by dist
end

if Limits > 0,
  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   double((DD(:,5)<100)), 'filled'); % color by true/near
end


%  scatter3(DD(:,10),DD(:,11),DD(:,12)-0.0001*DD(:,5), s,   DD(:,8), 'filled'); % color by angle
%  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   DD(:,7), 'filled'); % color by oxygen atom
%  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   abs(DD(:,12)), 'filled'); % color oxygen by vertical displacement

  for i = 1:length(DD(:,10)),
    if s(i) > 1,
%      plot3(DD(i,[10 13]), DD(i,[11 14]), DD(i,[12 15]), 'k');
    end
  end

  caxis([2.3 4.5]);                        % color axis for distances
if Limits > 0,
  caxis([0 1.5]);
end
%  caxis([80 180]);

  zPlotStandardBase(v,1,0);                % plot base at the origin
  title(['Phosphate interactions with ' L{v}]);
  rotate3d on
  grid off
  axis equal

    switch v
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
      case 2,                         % Base C
              h   = [10 11 12 13];
              m   = [ 7  8  6  6];
      case 3,                         % Base G
              h   = [12 13 15 16];
              m   = [ 4  7 11 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              m   = [ 8  4  7];
    end

 if Limits > 0,
  for hh = 1:length(h),
    vv = StandardLoc(h(hh),:,v) - StandardLoc(m(hh),:,v);
    vang = atan2(vv(2),vv(1));               % angle v makes with x axis
    xo = StandardLoc(h(hh),1,v);
    yo = StandardLoc(h(hh),2,v);

    a1 = 180-AL;
    th = (vang - a1*pi/180):0.1:(vang + a1*pi/180);
    x = xo + DL(m(hh))*cos(th);
    y = yo + DL(m(hh))*sin(th);
    x = [xo x xo];
    y = [yo y yo];
    plot(x,y,'k');

    a1 = 180-nAL;
    th = (vang - a1*pi/180):0.1:(vang + a1*pi/180);
    x = xo + nDL(m(hh))*cos(th);
    y = yo + nDL(m(hh))*sin(th);
    x = [xo x xo];
    y = [yo y yo];
    plot(x,y,'g');
  end
 end

  switch v,
    case 1,     text(8,0,'2BPh');
                text(5,8,'6BPh');
                text(-4.4,7.2,'7BPh');
                text(-3.5,-4,'0BPh');
                axis([-6 9 -5 9]);
    case 2,     text(5,6.1,'6BPh');
                text(-3.9,8,'7BPh');
                text(-5.7,6.6,'8BPh');
                text(-7.5,4.7,'9BPh');
                text(-4,-4,'0BPh');
                axis([-8 7 -5 9]);
    case 3,     text(7,-4,'1BPh');
                text(8.5,3,'3BPh');
                text(7.3,4.8,'4BPh');
                text(6,6,'5BPh');
                text(-4.5,-4,'0BPh');
                axis([-6 9 -5 9]);
    case 4,     text(6.2,3,'5BPh');
                text(-6.9,5.5,'9BPh');
                text(-2.9,-4,'0BPh');
                axis([-8 7 -5 9]);

  end
  view(2)

  if v == 2,
%    colorbar('east')
  end

end

  colormap('default')
  map = colormap;
%  map(1,:) = [0 0 0];
  colormap(flipud(map));


%set(gcf,'Renderer','OpenGL');     % fast rotation
%set(gcf,'Renderer','zbuffer')
set(gcf,'Renderer','painters');    % makes nice PDF files

saveas(gcf,['Phosphate Interactions\BaseOxygenPhosphorus.fig'],'fig');
saveas(gcf,['Phosphate Interactions\BaseOxygenPhosphorus.png'],'png');
saveas(gcf,['Phosphate Interactions\BaseOxygenPhosphorus.pdf'],'pdf');

%  saveas(gcf,['Phosphate Interactions\PhosphateInteractions_' L{v} '.fig'],'fig')

%  saveas(gcf,['Phosphate interactions' filesep 'Phosphate with ' L{v} '.png'],'png');

% ----------------- Display bond length and angle by type of massive base
figure(7)
clf


% temporary redefinition!

MT = {'Self interactions with C6/C8','Carbon with no overlap','G 3BPh/5BPh C 7BPh/9BPh (overlap)','Nitrogen with no overlap'};

ICode{1} = [4 9 14 17];             % carbon C6/C8 
ICode{2} = [1 16];                % ring carbons, but not C6/C8
ICode{3} = [6 8 11 13];                 % ring (imino) nitrogen
ICode{4} = [2 3 5 10 15];         % amino nitrogen

% correct definition:

MT = {'Self interactions with C6/C8','Carbon (excluding all self interactions)','Imino Nitrogen G(N1) U(N3)','Amino Nitrogen A(N6) C(N4) G(N2)'};

ICode{1} = [4 9 14 17];             % carbon C6/C8 
ICode{2} = [1 8 16];                % ring carbons, but not C6/C8
ICode{3} = [13 15];                 % ring (imino) nitrogen
ICode{4} = [2 3 5 6 10 11];         % amino nitrogen


for v = 1:4,
  subplot(2,2,v)

  r = [];
  for k = 1:length(ICode{v}),
    if v == 2,                         % exclude self interactions
      r = [r (find((mod(D(:,5),100) == ICode{v}(k)) .* (D(:,2) ~= D(:,3))))'];
    elseif v == 1,                     % keep only self interactions
      r = [r (find((mod(D(:,5),100) == ICode{v}(k)) .* (D(:,2) == D(:,3))))']; 
    else
      r = [r (find((mod(D(:,5),100) == ICode{v}(k))))'];  
    end
                                       % append matches to this massive atom
  end

  if v == 2,                           % add in non-self 0BP pairs
    for k = 1:length(ICode{2}),
      r = [r (find((mod(D(:,5),100) == ICode{2}(k)) .* (D(:,2) ~= D(:,3))))'];
    end
  end

%  if v == 1,                           % add in self C5-H5 pairs
  if v == 10,                           % add in self C5-H5 pairs
    for k = 1:length(ICode{1}),
      r = [r (find((mod(D(:,5),100) == ICode{1}(k)) .* (D(:,2) == D(:,3))))'];
    end
  end

  if (v == 1) && (length(r) > 4000),   % reduce number of self interactions
    y = rand(1,length(r));
    [w,i] = sort(y);
    r = r(i(1:4000));
  end

  DD = D(r,:);                         % use only these lines of data

%  s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
%  scatter3(DD(:,9),DD(:,8), DD(:,12), s,   DD(:,4), 'filled');

  s = 4*(DD(:,17) == 1) + 0*(DD(:,17) == 0);% best oxygen is large

%  scatter(DD(:,9),DD(:,8), s, abs(DD(:,12)), 'filled'); % color by vert displ
%  scatter(DD(:,9),DD(:,8), s,   DD(:,4), 'filled'); % color by base

  i = find((DD(:,7) == 1) .* (DD(:,17) == 1));
  scatter3(DD(i,9),DD(i,8), rand(length(i),1), 4,   'g', 'filled'); % color by oxygen atom
hold on
  i = find((DD(:,7) == 2) .* (DD(:,17) == 1));
  scatter3(DD(i,9),DD(i,8), rand(length(i),1), 4,   'b', 'filled'); % color by oxygen atom
  i = find((DD(:,7) >  2) .* (DD(:,17) == 1));
  scatter3(DD(i,9),DD(i,8), rand(length(i),1), 4,   'r', 'filled'); % color by oxygen atom

  view(2)
  grid off

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 4.5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 4.0;                           % near category

  hold on
  if v <=2,
    plot3([1 CarbonDist CarbonDist], [AL AL 180], [1 1 1], 'k', 'LineWidth', 1.5);
    plot3([1 nCarbonDist nCarbonDist], [nAL nAL 180], [1 1 1], 'k', 'LineWidth', 1.0);
  else
    plot3([1 NitrogenDist NitrogenDist], [AL AL 180], [1 1 1], 'k', 'LineWidth', 1.5);
    plot3([1 nNitrogenDist nNitrogenDist], [nAL nAL 180], [1 1 1], 'k', 'LineWidth', 1.0);
  end

  title([MT{v}]);
  if any(v==[1 3]),
    ylabel('Hydrogen bond angle');
  end
  if any(v==[3 4]),
    xlabel('Distance from hydrogen donor');
  end
  axis([2 5 90 180]);
  caxis([1 4]);
end

set(gcf,'Renderer','painters');    % makes nice PDF files

saveas(gcf,['Phosphate Interactions\BasePhosphateParameters.fig'],'fig');
saveas(gcf,['Phosphate Interactions\BasePhosphateParameters.png'],'png');
saveas(gcf,['Phosphate Interactions\BasePhosphateParameters.pdf'],'pdf');

return

% ----------------- Distance between centers of bases
figure(5)
clf
hist(D(:,16),30);
max(D(:,16));
title('Distance between centers of bases');

% ----------------- Display bond length and angle by nucleotide
figure(6)
clf
for v = 1:4,
  subplot(2,2,v)

  r = find(D(:,4) == v);               % select the current nucleotide code
  DD = D(r,:);                         % use only these lines of data

% s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
  s = 3*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

%  scatter(DD(:,9),DD(:,8), s,   DD(:,6), 'filled');  % color by massive atom
  scatter(DD(:,9),DD(:,8), s,   DD(:,7), 'filled');  % color by oxygen atom
  hold on
  title(['Oxygen location parameters for ', L{v}]);
  axis([2 5 100 180]);
  caxis([0 5]);
end



% -------------- Histogram of vertical displacement of phosphorus
figure(8)
clf
hist(D(:,15),30);
title('Vertical displacement of phosphorus for BP and near BP');

% -------------- Display angle and distance to phosphorus atom by type of massive atom
figure(9)
clf
MT = {'Carbon (excluding C6/C8 self interactions)','self interaction with C6/C8','Ring Nitrogen','Amino Nitrogen'};


ICode{1} = [1 8 16];                % ring carbon, not C6/C8
ICode{2} = [4 9 14 17];             % carbon C6/C8 self interactions 
ICode{3} = [13 15];                 % ring nitrogen
ICode{4} = [2 3 5 6 10 11];         % amino nitrogen

for v = 1:4,
  subplot(2,2,v)

  r = [];
  for k = 1:length(ICode{v}),
    if v == 2                          % keep only self interactions
      r = [r (find((mod(D(:,5),100) == ICode{v}(k)) .* (D(:,2) == D(:,3))))'];  
    else
      r = [r (find(mod(D(:,5),100) == ICode{v}(k)))'];  
    end
                                       % append matches to this massive atom
  end

  if v == 1,                           % add in non-self 0BP pairs
    for k = 1:length(ICode{2}),
      r = [r (find((mod(D(:,5),100) == ICode{2}(k)) .* (D(:,2) ~= D(:,3))))'];
    end
  end

  DD = D(r,:);                         % use only these lines of data

%  s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
  s = 4*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

  scatter(DD(:,20),DD(:,19), s,   DD(:,4), 'filled');
  hold on
%  plot([1 3.6 3.6], [150 150 180], 'k')
%  plot([1 4.0 4.0], [120 120 180], 'k')
  title(['Phosphorus location parameters for ', MT{v}]);
  axis([3 6 100 180]);
  caxis([0 5]);
end

% return

u = unique(D(:,[1 2 3 6]),'rows');
c = zeros(5,10);
CarbonDist    = 4.0;                           % max massive - oxygen distance
NitrogenDist  = 3.5;                           % max massive - oxygen distance
DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;
for h = 1:length(u(:,1)),
  i = find( (D(:,1)==u(h,1)) .* (D(:,2)==u(h,2)) .* (D(:,3)==u(h,3)) .* (D(:,6)==u(h,4)));
    switch D(i(1),4),
      case 1,                         % Base A
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
      case 2,                         % Base C
              m   = [ 7  8  6  6];
      case 3,                         % Base G
              m   = [ 4  7 11 11];
      case 4,                         % Base U
              m   = [ 8  4  7];
    end
  a = sum(D(i,5) < 100);
  b = sum((D(i,5) >= 100) .* (D(i,9) < DL(m(u(h,4)))));   % near but OK dist
  d = sum((D(i,5) >= 100) .* (D(i,9) >= DL(m(u(h,4)))));  % near b/c large dist
  c(a+1,b+1) = c(a+1,b+1) + 1;
  c(a+1,5+d+1) = c(a+1,5+d+1) + 1;

  if b == 2,
    for y = 1:length(i),
            f = D(i(y),1);
            N1 = File(f).NT(D(i(y),2));
            N2 = File(f).NT(D(i(y),3));

            fprintf('%6s base %s%5s BPcode %3d %4s phosphate donor %s%5s length %6.2f angle %6.2f\n', File(f).Filename, N1.Base, N1.Number, D(i(y),5), zBasePhosphateText(D(i(y),5)), N2.Base, N2.Number, D(i(y),9), D(i(y),8));
    end

    j = find( (D(:,1)==u(h,1)) .* (D(:,2)==u(h,2)) .* (D(:,3)==u(h,3)) .* (D(:,6)~=u(h,4)));
    if length(j) > 0,
      fprintf('This pair also makes the following interaction:\n');
      for y = 1:length(j),
            f = D(j(y),1);
            N1 = File(f).NT(D(j(y),2));
            N2 = File(f).NT(D(j(y),3));

            fprintf('%6s base %s%5s BPcode %3d %4s phosphate donor %s%5s length %6.2f angle %6.2f\n', File(f).Filename, N1.Base, N1.Number, D(j(y),5), zBasePhosphateText(D(j(y),5)), N2.Base, N2.Number, D(j(y),9), D(j(y),8));
      end
    end  
    fprintf('\n');

  end

end

fprintf('                                   Near inter w/ OK distance.   Near inter. w/ too large dist.\n');
fprintf('                                     ');
for b = 0:4,
  fprintf(' %4d', b);
end
for b = 0:4,
  fprintf(' %4d', b);
end
fprintf(' Total\n');
for a = 0:4,
  fprintf('%4d oxygens making true interaction ', a);
  for b = 0:9,
    fprintf(' %4d', c(a+1,b+1));
  end
  fprintf(' %4d\n', sum(c(a+1,:))/2);
end
fprintf('Total                                ');
for b = 0:9,
  fprintf(' %4d', sum(c(:,b+1)));
end
fprintf(' %4d\n', sum(sum(c))/2);

return

r = find((D(:,10) < 3) .* (D(:,10) > 0) .* (D(:,11) < 3) .* (D(:,11) > 0));

[a,b,c] = unique(D(r,[1 2 3]),'rows');
Search = [D(r(b),2) D(r(b),3) D(r(b),1)];
xDisplayCandidates(File,Search)
