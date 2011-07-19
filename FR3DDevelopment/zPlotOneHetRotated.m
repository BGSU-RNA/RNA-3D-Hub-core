% zPlotOneHetRotated(Het,ViewParam,R,S) rotates Het and displays with ViewParam

function [void] = zPlotOneHetRotated(Het,ViewParam,R,S)

L           = length(Het.Loc(:,1));         % Number of Het atoms
NewHet.Loc   = (Het.Loc   - ones(L,1) *S) * R; % rotated into position
NewHet.Unit  = Het.Unit;
NewHet.Number= Het.Number;

if isfield(Het,'Beta'),
  NewHet.Beta  = Het.Beta;
end

zPlotOneHet(NewHet,ViewParam)
