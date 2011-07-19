% zCoplanarStats analyzes basepaired nucleotides in File to determine observed values of the angle between the planes of the bases (Pair.PlaneAng), the gap between their interacting edges (Pair.Gap), and the minimum distance between the bases (Pair.MinDist).  


% File = zAddNTData('1s72');

Stats = [];

for f = 1:length(File),
  E = File(f).Edge;
  [i,j] = find( (abs(E) > 0) .* (abs(E) < 14) );
  for k = 1:length(i),
    N1 = File(f).NT(i(k));
    N2 = File(f).NT(j(k));


%  Pair.PlaneAng = acos(abs(N1.Rot(:,3)'*N2.Rot(:,3)))*57.29577951308232; 
                                             % angle between planes

  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  d = zDistance(N2.Fit(1:Lim(2,N2.Code),:), N1.Center); 
                                           % distances to base 1 center
  [y,m] = min(d);                          % identify the closest atom
  m = m(1);                                % in case of a tie, use the first
  Pair.Gap = N1.Rot(:,3)'*(N2.Fit(m,:)-N1.Center)';% height above plane of 1

  Pair.MinDist = min(min(zDistance(N1.Fit,N2.Fit)));

  % ------------------------ check for coplanarity

  Pair.Coplanar = 0;                      % default value, not coplanar

  % Criteria for possibly being coplanar:
  %   Pair.Gap must be < 2 Angstroms
  %   Pair.MinDist must be < 4.5 Angstroms
  %   Angle between center-center vector and each normal must be > 60 degrees
  %   Angle between normal vectors must be < 45 degrees


    v  = N1.Center - N2.Center;           % vector from center to center
    v  = v / norm(v);                     % normalize

    dot1 = abs(v * N1.Rot(:,3));          % to calculate angle: v and normal
    dot2 = abs(v * N2.Rot(:,3));
    dot3 = abs(N1.Rot(:,3)' * N2.Rot(:,3));

    yy = 0.5;                             % cos(60) = 0.5
    yyy = 1/sqrt(2);                      % cos(45) = 1/sqrt(2)

      d = zDistance(N1.Fit(1:Lim(2,N1.Code),:), N2.Center); 
                                           % distances to base 2 center
      [y,m] = min(d);                      % identify the closest atom
      m = m(1);                            % in case of a tie, use the first
      Gap2 = N2.Rot(:,3)'*(N1.Fit(m,:)-N2.Center)';% height above plane of 1

        Pair.Coplanar = min([(2-abs(Pair.Gap))/2 (2-abs(Gap2))/2 (yy-dot1)/yy (yy-dot2)/yy (yyy-dot3)/yyy min(1,4.5-Pair.MinDist)]);


%      Stats = [Stats; [(2-abs(Pair.Gap))/2 (2-abs(Gap2))/2 (yy-dot1)/yy (yy-dot2)/yy (dot3-yyy)/yyy min(1,4.5-Pair.MinDist)]];

      Stats = [Stats; [abs(Pair.Gap) abs(Gap2) dot1 dot2 dot3 Pair.MinDist]];

  end
end

names{1} = 'abs(Pair.Gap)';
names{2} = 'abs(Gap2)';
names{3} = 'dot1';
names{4} = 'dot2';
names{5} = 'dot3';
names{6} = 'Pair.MinDist';




for v = 1:6,
  figure(v)
  clf
  if v == 6,
    i = find(Stats(:,v)<4);
  else
    i = 1:length(Stats(:,1));
  end
  N = length(i);
  hist(Stats(i,v),30)
  title(names{v});
  hold on
  plot(quantile(Stats(i,v),0.10),0,'r*');
  plot(quantile(Stats(i,v),0.90),0,'r*');
  plot(quantile(Stats(i,v),0.97),0,'c*');
  plot(quantile(Stats(i,v),0.03),0,'c*');
  plot(quantile(Stats(i,v),0.96),0,'c*');
  plot(quantile(Stats(i,v),0.04),0,'c*');
  plot(quantile(Stats(i,v),0.50),0,'m*');
  plot(quantile(Stats(i,v),0.70),0,'g*');
  plot(quantile(Stats(i,v),0.30),0,'g*');

  if v ~= 5,
    fprintf('Contribution to Coplanar: 1+(%s-%7.4f)*(%7.4f)\n', names{v}, quantile(Stats(i,v),0.70), 0.5/(quantile(Stats(i,v),0.70)-quantile(Stats(i,v),0.90)));

    fprintf('Upper limit to be ncp:    %s <= %7.4f\n', names{v}, (-quantile(Stats(i,v),0.70)+2*quantile(Stats(i,v),0.90)));

    L = (-quantile(Stats(i,v),0.70)+2*quantile(Stats(i,v),0.90));
    fprintf('Percentage below this:    %7.4f\n', 100*sum(Stats(i,v) < L)/N);
    fprintf('97th percentile:          %7.4f\n', quantile(Stats(i,v),0.97)); 

  else
    fprintf('Contribution to Coplanar: 1+(%s-%7.4f)*(%7.4f)\n', names{v}, quantile(Stats(i,v),0.30), 0.5/(quantile(Stats(i,v),0.30)-quantile(Stats(i,v),0.10)));
    fprintf('Lower limit to be ncp:    %s >= %7.4f\n', names{v}, (-quantile(Stats(i,v),0.30)+2*quantile(Stats(i,v),0.10)));
    L = (-quantile(Stats(i,v),0.30)+2*quantile(Stats(i,v),0.10));
    fprintf('Percentage above this:    %7.4f\n', 100*sum(Stats(i,v) < L)/N);
    fprintf('3rd percentile:           %7.4f\n', quantile(Stats(i,v),0.03)); 

  end

  fprintf('\n');

  q(v) = quantile(Stats(i,v),0.97);
  p(v) = quantile(Stats(i,v),0.90);

end

q(5) = quantile(Stats(i,5),0.03);
B97 = (Stats(:,1) <= q(1)) .* (Stats(:,2) <= q(2)) .* (Stats(:,3) <= q(3)) .* (Stats(:,4) <= q(4)) .* (Stats(:,5) >= q(5)) .* (Stats(:,6) <= q(6)) ;

sum(B97)/N  

p(5) = quantile(Stats(i,5),0.10);
B97 = (Stats(:,1) <= p(1)) .* (Stats(:,2) <= p(2)) .* (Stats(:,3) <= p(3)) .* (Stats(:,4) <= p(4)) .* (Stats(:,5) >= p(5)) .* (Stats(:,6) <= p(6)) ;

sum(B97)/N  


names{5} = '-dot3';

VarNames{1} = 'Gap1Val';
VarNames{2} = 'Gap2Val';
VarNames{3} = 'dot1Val';
VarNames{4} = 'dot2Val';
VarNames{5} = 'dot3Val';
VarNames{6} = 'MinDistVal';

fprintf('       %% Code generated by zCoplanarStats.m\n');

for v = 1:6,
  if v == 6,
    i = find(Stats(:,v)<4);
  else
    i = 1:length(Stats(:,1));
  end
  N = length(i);

  if v == 5,
    q70 = quantile(-Stats(i,v),0.70);
    q90 = quantile(-Stats(i,v),0.90);
    q97 = quantile(-Stats(i,v),0.97);
  else
    q70 = quantile(Stats(i,v),0.70);
    q90 = quantile(Stats(i,v),0.90);
    q97 = quantile(Stats(i,v),0.97);
  end

  fprintf('      if %s < %7.4f,               %% 70th percentile\n', names{v}, q70);
  fprintf('        %s = 1;\n', VarNames{v});
  fprintf('      elseif %s < %7.4f,           %% 90th percentile\n', names{v}, q90);
  fprintf('        %s = 1+(%s-%7.4f)*(%7.4f);\n', VarNames{v}, names{v}, q70, 0.5/(q70-q90));
  fprintf('      elseif %s < %7.4f,           %% 97th percentile\n', names{v}, q97);
  fprintf('        %s = 0.5+(%s-%7.4f)*(%7.4f);\n', VarNames{v}, names{v}, q90, 0.5/(q90-q97));
  fprintf('      else\n');
  fprintf('        %s = 0;\n', VarNames{v});
  fprintf('      end\n');
  fprintf('\n');
end

break

for dot3 = 0.7:0.01:1,
    if -dot3 < -0.9509,
      dot3Val = 1;
    elseif -dot3 < -0.8835,
      dot3Val = 1+(-dot3--0.9509)*(-7.4217);
    elseif -dot3 < -0.7757,
      dot3Val = 0.5+(-dot3--0.8835)*(-4.6390);
    else
      dot3Val = 0;
    end
[dot3 dot3Val]
end
