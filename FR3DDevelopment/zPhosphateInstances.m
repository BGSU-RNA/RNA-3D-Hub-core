
%F = zAddNTData('2avy');

f = 2;

F = File(f);

[i,j,k] = find(F.BasePhosphate);


 for c = 1:17,
  g = find(k == c);

%for c = [1:7 18 8:12 19 13:17],
%  g = find((k == c) .* (i ~= j));          % non-self interactions only

  [y,z]= sort(i(g));

  g = g(z);

  for a = 1:min(5,length(g)),
%  for a = 1:min(500,length(g)),
    N1 = F.NT(i(g(a)));
    N2 = F.NT(j(g(a)));

    J = N1.Hierarchy{3};
    if ~isempty(strfind(J,'HL')),
      H = 'Hairpin    ';
    elseif ~isempty(strfind(J,'H')),
      H = 'Helix      ';
    elseif ~isempty(strfind(J,'IL')),
      H = 'Internal   ';
    elseif ~isempty(strfind(J,'BL')),
      H = 'Bulge      ';
    elseif ~isempty(strfind(J,'J')),
      H = 'Junction   ';
    elseif ~isempty(strfind(J,'3_p')),
      H = '3 prime end';
    elseif ~isempty(strfind(J,'5_p')),
      H = '5 prime end';
    end

    H1 = H;

    J = N2.Hierarchy{3};
    if ~isempty(strfind(J,'HL')),
      H = 'Hairpin    ';
    elseif ~isempty(strfind(J,'H')),
      H = 'Helix      ';
    elseif ~isempty(strfind(J,'IL')),
      H = 'Internal   ';
    elseif ~isempty(strfind(J,'BL')),
      H = 'Bulge      ';
    elseif ~isempty(strfind(J,'J')),
      H = 'Junction   ';
    elseif ~isempty(strfind(J,'3_p')),
      H = '3 prime end';
    elseif ~isempty(strfind(J,'5_p')),
      H = '5 prime end';
    end

    H2 = H;

    if File(f).Crossing(i(g(a)),j(g(a))) < 2,   % number of cWW's crossed
      R = 'Local';
    else
      R = 'Long-range';
    end

    fprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s - %s - %s\n', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBasePhosphateText(F.BasePhosphate(i(g(a)),j(g(a))),0), H1, H2, R);
  end
end

