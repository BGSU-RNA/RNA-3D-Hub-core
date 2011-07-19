% zOtherPhosphateInteractions calculates possible base-phosphate interactions not mediated by a base hydrogen

    switch N1.Code
      case 1,                         % Base A
              m   = [ 4  3  10];      % rows of the corresponding massive atoms
      case 2,                         % Base C
              m   = [ 7  8  6  6];
      case 3,                         % Base G
              m   = [ 4  7 11 11];
      case 4,                         % Base U
              m   = [ 8  4  7];
    end

    dis = zDistance(N1.Fit(m,:), N2.Sugar(p,:)); % distances between mass & O's
    nearDL = nDL(m)' * ones(1,4);     % limits to compare to
    dis = dis .* (dis < nearDL);      % massive-oxygen pairs close enough

    g = [];                           % internal classification number
    Angle = [];
    Dist  = [];

    for mm = 1:length(m),             % massive atom to consider
     ppp = find(dis(mm,:));            % oxygens close enough to consider
     for n = 1:length(ppp),            % loop through potential oxygens
      Angle(n)=zAngle(N1.Fit(m(mm),:),N1.Fit(h(mm),:),N2.Sugar(p(ppp(n)),:));
                                      % base massive - hydrogen - oxygen angle
      Dist(n) = dis(mm,ppp(n));        % distance
     end

     PAngle = zAngle(N1.Fit(m(mm),:),N1.Fit(h(mm),:),N2.Sugar(10,:));
                                  % base massive - hydrogen - phosphorus angle
     PDist  = zDistance(N1.Fit(m(mm),:),N2.Sugar(10,:));        % distance

     [u,v] = min(-Angle+60*Dist);     % order by quality of potential bond

     for n = 1:length(ppp),            % loop through potential oxygens
      if Angle(n) > nAL,              % good enough to be "near" base-phosph

        if ((Angle(n) > AL) && (Dist(n) < DL(m(mm)))) % true BP pair
          g = [g e(mm)];              % assign a non-near class.
          T = [T; [f i(k) j(k) e(mm) ppp(n)]];
        else
          g = e(mm) + 100;            % > 100 means "near"
        end

        if Verbose > 1,
          % store information for later display
          ox = (N2.Sugar(p(ppp(n)),:)-N1.Fit(1,:)) * N1.Rot; % oxygen displ

          a = [f i(k) j(k) N1.Code g(end) mm ppp(n) Angle(n) Dist(n) ox ph File(f).Distance(i(k),j(k)) (v(1)==n) -u(1) PAngle PDist];

% [ g(end) v(1) == n]

          % Columns:
          % 1  file number
          % 2  index of base
          % 3  index of nucleotide using phosphate
          % 4  code of base
          % 5  classification number for this massive-oxygen pair
          % 6  which massive atom is interacting
          % 7  which oxygen is interacting
          % 8  angle of interaction, in degrees
          % 9  distance from massive atom to oxygen, in Angstroms
          %10  displacement of oxygen atom relative to C1' of base
          %13  displacement of phophorus atom relative to C1' of base
          %16  distance between centers of the two bases
          %17  1 if this is the best oxygen for this hydrogen, 0 otherwise
          %18  approximate quality of the hydrogen bond
          %19  angle made by massive, hydrogen, phosphorus
          %20  distance from massive to phosphorus

          D = [D; a];                  % append data to data matrix

          if Verbose > 3,
            fprintf('%6s base %s%5s %3s BPcode %3d %4s phosphate donor %s%5s %13s length %6.2f angle %6.2f interaction %s', File(f).Filename, N1.Base, N1.Number, AtomNames{h(mm),N1.Code}, g(end), zBasePhosphateText(g(end)), N2.Base, N2.Number, Sugar{p(ppp(n))}, dis(mm,ppp(n)), Angle(n), zEdgeText(File(f).Edge(i(k),j(k))));
            if a(17) == 1,
              fprintf('Best\n');
            else
              fprintf('\n');
            end
          end
        end

      end

     end  % loop over potential oxygens
    end   % loop over massive atoms
