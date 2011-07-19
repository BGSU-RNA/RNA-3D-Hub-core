function [void] = zPlotOneHet(Het,ViewParam)

if isfield(ViewParam,'LineStyle'),
  LS = ViewParam.LineStyle;
else
  LS = '-';
end

if isfield(ViewParam,'LineThickness'),   
  LT = ViewParam.LineThickness;
else
  LT = 2.0;                              % doesn't work for some reason!
end

if isfield(ViewParam,'LabelBases'),
  LB = ViewParam.LabelBases;
else
  LB = 0;
end

if isfield(ViewParam,'LabelSugar'),
  LSugar = ViewParam.LabelSugar;
else
  LSugar = 0;
end

if isfield(ViewParam,'ShowBeta'),
  ShowBeta = ViewParam.ShowBeta;
else
  ShowBeta = 0;
end

gray = 0.5*[1 1 1];

col = [1 0 1];

if isfield(ViewParam,'Color'),
  if length(ViewParam.Color) == 3,
    col = ViewParam.Color;
  end
end

if strcmp(LS,'-.'),
  col = 0.7*col;
  gray = 0.7*gray;
  LS = '-';
end

bc = gray;

hold on 

X  = Het.Loc;

scatter3(X(:,1),X(:,2),X(:,3),10,col,'filled');

if LB > 0,
  text(X(1,1)+0.1,X(1,2),X(1,3)+0.1,[Het.Unit Het.Number],'fontweight','bold','FontSize',LB);
%  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[Het.Unit Het.Number],'fontweight','bold','FontSize',LB);
end

if ShowBeta > 0 && isfield(Het,'Beta'),
  for j = 1:min(12,length(Het.Beta(:,1))),
    if Het.Beta(j,1) < Inf,
      text(Het.Sugar(j,1),Het.Sugar(j,2),Het.Sugar(j,3),num2str(round(Het.Beta(j,1))));
    end
  end
  for j = 13:length(Het.Beta),
    k = j - 12;
    if Het.Beta(k,1) < Inf,
      text(Het.Loc(j-12,1),Het.Loc(j-12,2),Het.Loc(j-12,3),num2str(round(Het.Beta(j,1))));
    end
  end
end
      

