% zPhosCompare compares the base-phosphate interaction parameters for
% equivalent 3D structures.

conv = [1 2 3 3];
D(:,7) = conv(D(:,7));                 % O1P and O2P are equivalent

L = {'A','C','G','U'};

% display interacting oxygens and phosphorus atoms together with bases


k = find( (D(:,21) ~= D(:,22)) .* (D(:,17) == 1));
DD = D(k,:);

figure(1)

Lett = 'ACGU';

m  = 2;
b  = 2;

i = find( (DD(:,6) == m) .* (DD(:,4) == b) );
DD = DD(i,:);

while 0 < 1,
  clf

  j = ceil(rand*length(DD(:,1)));               % choose a row at random
  i = find( (DD(:,21) == DD(j,21)) .* (DD(:,22) == DD(j,22)));  % same nucls

  DD(i,:)

  scatter(DD(i,9),DD(i,8), 'filled');  % color by oxygen atom
  hold on
  title(['Oxygen location parameters for ' Lett(DD(j,4)) num2str(DD(j,21))]);
  axis([2 5 100 180]);
  caxis([0 5]);

  drawnow
  pause
end




break

figure(1);
clf

for v = 1:4,
%  figure(v)
%  clf
  subplot(2,2,v);

%  r = find(D(:,4) == v);               % select the current nucleotide code

  r = find((D(:,4) == v) .* (abs(D(:,12)) < 3)); % select nucl and vert displ

  DD = D(r,:);                         % use only these lines of data

  s = 10*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
%  s = 10*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

  scatter3(DD(:,13),DD(:,14),DD(:,15), s, 0*DD(:,9), 'filled');
  hold on

%  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   DD(:,9), 'filled'); % color by dist
%  scatter3(DD(:,10),DD(:,11),DD(:,12)-0.0001*DD(:,5), s,   DD(:,8), 'filled'); % color by angle
%  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   DD(:,7), 'filled'); % color by oxygen atom
  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   abs(DD(:,12)), 'filled'); % color by vertical displacement


  for i = 1:length(DD(:,10)),
    if s(i) > 1,
%      plot3(DD(i,[10 13]), DD(i,[11 14]), DD(i,[12 15]), 'k');
    end
  end

  switch v,
    case 1,     text(8,0,'2BP');
                text(5,8,'6BP');
                text(-3,7.2,'7BP');
                text(-3.5,-3,'0BP');
    case 2,     text(5,6.1,'6BP');
                text(-2.3,8,'7BP');
                text(-4.5,6.3,'8BP');
                text(-5.8,4.7,'9BP');
                text(-4,-2.8,'0BP');
    case 3,     text(8,-2.4,'1BP');
                text(8.5,3,'3BP');
                text(7.3,4.8,'4BP');
                text(6,6,'5BP');
                text(-4.5,-2.9,'0BP');
    case 4,     text(6.2,3,'5BP');
                text(-3.9,6,'9BP');
                text(-2.9,-2,'0BP');
  end

  map = colormap;
  map(1,:) = [0 0 0];
  colormap(map);

  caxis([0 4]);
%  caxis([80 180]);

  zPlotStandardBase(v,1,0);                % plot base at the origin
  title(['Phosphate interactions with ' L{v}]);
  rotate3d on
  grid off
  axis([-6 9 -5 9]);
%  axis equal
  view(2)
%  saveas(gcf,['Phosphate Interactions\PhosphateInteractions_' L{v} '.fig'],'fig')

%  saveas(gcf,['Phosphate interactions' filesep 'Phosphate with ' L{v} '.png'],'png');

set(gcf,'Renderer','OpenGL');     % fast rotation
%set(gcf,'Renderer','zbuffer')
%set(gcf,'Renderer','painters');    % makes nice PDF files

end

figure(5)
clf
hist(D(:,16),30);
max(D(:,16));

% display bond length and angle by nucleotide

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

% display bond length and angle by type of massive base

figure(7)
clf

MT = {'Carbon (excluding all self interactions)','all self interaction with C5/C6/C8','Ring Nitrogen','Amino Nitrogen'};

ICode{1} = [1 8 16];                % ring carbons, but not C6/C8
ICode{2} = [4 9 14 17];             % carbon C6/C8 
ICode{3} = [13 15];                 % ring nitrogen
ICode{4} = [2 3 5 6 10 11];         % amino nitrogen

for v = 1:4,
  subplot(2,2,v)

  r = [];
  for k = 1:length(ICode{v}),
    if v == 1,                         % exclude self interactions
      r = [r (find((mod(D(:,5),100) == ICode{v}(k)) .* (D(:,2) ~= D(:,3))))'];
    elseif v == 2,                     % keep only self interactions
      r = [r (find((mod(D(:,5),100) == ICode{v}(k)) .* (D(:,2) == D(:,3))))']; 
    else
      r = [r (find((mod(D(:,5),100) == ICode{v}(k))))'];  
    end
                                       % append matches to this massive atom
  end

  if v == 1,                           % add in non-self 0BP pairs
    for k = 1:length(ICode{2}),
      r = [r (find((mod(D(:,5),100) == ICode{2}(k)) .* (D(:,2) ~= D(:,3))))'];
    end
  end

  if v == 2,                           % add in self C5-H5 pairs
    for k = 1:length(ICode{1}),
      r = [r (find((mod(D(:,5),100) == ICode{1}(k)) .* (D(:,2) == D(:,3))))'];
    end
  end

  if (v == 2) && (length(r) > 5000),
    y = rand(1,length(r));
    [w,i] = sort(y);
    r = r(i(1:5000));
  end

  DD = D(r,:);                         % use only these lines of data

%  s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
  s = 4*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

%  scatter3(DD(:,9),DD(:,8), DD(:,12), s,   DD(:,4), 'filled');

%  scatter(DD(:,9),DD(:,8), s, abs(DD(:,12)), 'filled'); % color by vert displ
%  scatter(DD(:,9),DD(:,8), s,   DD(:,4), 'filled'); % color by base
  scatter(DD(:,9),DD(:,8), s,   DD(:,7), 'filled'); % color by oxygen atom

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 4.5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 4.0;                           % near category

  hold on
  if v <=2,
    plot([1 CarbonDist CarbonDist], [AL AL 180], 'k');
    plot([1 nCarbonDist nCarbonDist], [nAL nAL 180], 'k');
  else
    plot([1 NitrogenDist NitrogenDist], [AL AL 180], 'k');
    plot([1 nNitrogenDist nNitrogenDist], [nAL nAL 180], 'k');
  end

  title(['Oxygen location parameters for ', MT{v}]);
  axis([2 4.6 105 180]);
  caxis([1 4]);
end

figure(8)
clf
hist(D(:,15),30);
title('Vertical displacement of phosphorus for BP and near BP');

% display angle and distance to phosphorus atom by type of massive atom

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

break

r = find((D(:,10) < 3) .* (D(:,10) > 0) .* (D(:,11) < 3) .* (D(:,11) > 0));

[a,b,c] = unique(D(r,[1 2 3]),'rows');
Search = [D(r(b),2) D(r(b),3) D(r(b),1)];
xDisplayCandidates(File,Search)
