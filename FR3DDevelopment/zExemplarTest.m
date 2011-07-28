
% load PairExemplars

[s,t] = size(Exemplar);

VP.Sugar = 1;
VP.AtOrigin = 1;

OldClass = 0;

for j = 1:t,
  for i = s:-1:1,
    if ~isempty(Exemplar(i,j).NT1),
      pc = Exemplar(i,j).NT1.Code + 4*(Exemplar(i,j).NT2.Code-1);
      if j ~= pc,
        disp('The following exemplar has the wrong bases!');
        [i j pc]
        Exemplar(i,j)
      end

      F.NT(1) = Exemplar(i,j).NT1;
      F.NT(2) = Exemplar(i,j).NT2;

      Exemplar(i,j)

warning off

      figure(1)

      if fix(Exemplar(i,j).Class) ~= fix(OldClass),
        reset(gcf)
        clf
        cla
      end
      zDisplayNT(F,1:2,VP);
      title(['Class ' num2str(Exemplar(i,j).Class) ' ' zEdgeText(Exemplar(i,j).Class) ' ' Exemplar(i,j).NT1.Base Exemplar(i,j).NT2.Base]);
      view(2)

      figure(2)
      if fix(Exemplar(i,j).Class) ~= fix(OldClass),
        reset(gcf)
        clf
        cla
      end
      zDisplayNT(F,1:2,VP);
      title(['Class ' num2str(Exemplar(i,j).Class) ' ' zEdgeText(Exemplar(i,j).Class) ' ' Exemplar(i,j).NT1.Base Exemplar(i,j).NT2.Base]);
      view([48 0])
      pause

      OldClass = Exemplar(i,j).Class;

    end
  end
end
