% zPlotStandardBase(code,textoption) plots base with code,
% where A=1, C=2, G=3, U=4 and textoption:
% textoption = 0 - no text
% textoption = 1 - label atoms

% centeroption controls where the base is centered
% centeroption 0 - glycosidic atom at the origin (default)
% centeroption 1 - geometric center at the origin

% function [void] = zPlotBasePhosphateHydrogens(code,textoption,centeroption)

figure(1)
clf

Lett = 'ACGU';

shifts = zeros(19,2);                  % xy shifts for labels
shifts(8,:) = [0 0.1];
shifts(16,:) = [0 -0.1];
shifts(14,:) = [0 0.1];
shifts(4,:) = [0 -0.1];
shifts(9,:) = [0 0.1];
shifts(17,:) = [0 -0.1];

for code = 4:-1:1,

centeroption = 0;

textoption = 0;

  zStandardBases

  VP.Sugar = 0;

  BaseNames = 'ACGU';

  L = Lim(2,code);
  Q = StandardLoc(1:L,:,code);

  if centeroption == 1,
    M = Lim(1,code);
    Q = Q - ones(L,1)*mean(Q(1:M,:));
  end

  NT.Code = code;
  NT.Fit = Q;

switch code
  case 1,
    col = [1 0 0];   % A is red 
  case 2,
    col = [1 0.8 0]; % C is yellow
  case 3,
    col = [0 1 0];   % G is green
  case 4, 
    col = [0 0 1];   % U is blue
  otherwise,
    col = [0 0 0];
end

  VP.Color = col;                 % dim the bases

  zPlotOneNT(NT,VP);

    switch code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

  for z = 1:length(h),
    dx = Q(h(z),1)-Q(m(z),1);
    dy = Q(h(z),2)-Q(m(z),2);

    d = [dx dy];
    d = d / norm(d);

    if AtomNames{m(z),code}(1) == 'N',
      r = 2.9;
    elseif AtomNames{m(z),code}(1) == 'C',
      r = 3.3;
    else
      fprintf('Unrecognized heavy atom\n');
      code 
      m(z)
      AtomNames{m(z),code}
      r = 6;
    end

    ox = Q(m(z),1:2) + r*d;

    plot([Q(h(z),1) ox(1)],[Q(h(z),2) ox(2)], ':k', 'LineWidth', 2);

    plot(ox(1),ox(2), '.k', 'MarkerSize', 14);

%    plot([Q(h(z),1) Q(m(z),1)],[Q(h(z),2) Q(m(z),2)], 'Color', col, 'LineWidth', 2);
    plot(Q(h(z),1),Q(h(z),2), '.', 'Color', col, 'MarkerSize', 40);

    Tex = [' ' Lett(code) ' ' zBasePhosphateText(e(z))];

    switch e(z),
    case {11, 13}
      Tex = [Tex ', 4BPh'];
    case {6, 8},
      Tex = [Tex ', 8BPh'];
    end

    text(ox(1)+shifts(e(z),1),ox(2)+shifts(e(z),2),Tex);
%    text(ox(1)+shifts(e(z),1),ox(2)+shifts(e(z),2),[' ' num2str(e(z)) ' ' zBasePhosphateText(e(z))]);

  end


  H = [13 9 14 10];

  Z = [Q(H(code),:); Q(1,:)];
  k = [1 2]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',0.5*[1 1 1],'LineWidth',2,'LineStyle','-');
  if textoption > 0,
    hold on

    for j=1:L,
      text(Q(j,1), Q(j,2), 0.1, AtomNames{j,code},'FontSize',6);
    end
    title(['Standard ' BaseNames(code)])
  end

  axis off
  axis([-3.5    7.7   -1.8819    6.6345])
  axis equal

end

saveas(gcf,'BasePhosphateHydrogens.png','png');
