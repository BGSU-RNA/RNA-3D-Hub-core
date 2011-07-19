
figure(1)
clf
for v = 1:4,
  subplot(2,2,v)
  zPlotStandardBase(v,1)
  axis equal
end

break

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
