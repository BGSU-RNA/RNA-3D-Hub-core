
load PairExemplars

n1 = 3;                       % code of primary nucleotide
n2 = 2;                       % code of secondary nucleotide

if n1 == 5,
  A = [1 2 5 6   -3 -4 7 8 -5 -6 -9 11 -11 12];
else
  A = [1 2 3 4 5 6 -3 -4 7 8 9 10 -5 -6 -9 -10 11 -11 12 -12];
end

VP.AtOrigin = 1;
VP.Sugar = 1;
VP.LabelBases = 0;
VP.LineThickness = 3;

RNA = 'ACGU';

for i = 1:length(A),
  [NT1,NT2,E] = zGetExemplar(A(i),n1,n2,Exemplar);
  if isempty(E.Filename),
    A(i)
  else
    F.NT(1) = NT1;
    F.NT(2) = NT2;
    figure(1)
    clf
    zDisplayNT(F,'all',VP);
    view(2)
    grid off
    axis([-10 10 -10 10]);
    axis off
    n = zEdgeText(A(i));
    n = [n(1:2) ' / ' n(3)];
    n = strrep(n,'c','cis ');
    n = strrep(n,'t','trans ');
    n = strrep(n,'W','Watson-Crick ');
    n = strrep(n,'H','Hoogsteen ');
    n = strrep(n,'S','Sugar edge ');
    n = [RNA(n1) RNA(n2) ' ' n];
    title(n,'fontsize',20);
    j = num2str(i);
    if i < 10,
      j = ['0' j];
    end
    saveas(gca,[j '_' RNA(n1) RNA(n2) '_' zEdgeText(A(i)) '.png']);
    pause
  end
end


