% zExemplarIDI(Class,Code1,Code2) returns the 4x4 matrix of IDI values for pair Code1,Code2 in geometric family Class.  You can use text for these, for example, zExemplarIDI('cWW','A','U');  For speed, load PairExemplars and pass in ExemplarIDI.

function [IDI] = zExemplarIDI(Class,Code1,Code2,ExemplarIDI)

if nargin < 4,
  load PairExemplars
end

if strcmp(class(Class),'char'),
  Class = xGetEdgeNums(Class);
  Class = Class(1);
end

if strcmp(class(Code1),'char'),
  Code1 = pLtoN(Code1);
end

if strcmp(class(Code2),'char'),
  Code2 = pLtoN(Code2);
end

Class = fix(Class);

if Class < 0,
  Class = Class + 100;
end

IDI = ExemplarIDI{Class,Code1,Code2};

return


% -------------------------------- Display exemplars from each family

VP.AtOrigin = 1;

for Class = 1:12,
  figure(Class)
  clf
  for Code1 = 1:4,
    for Code2 = 1:4,
      [NT1,NT2,E] = zGetExemplar(Class,Code1,Code2);
      if ~isempty(NT1.Code),
        subplot(4,4,4*(Code1-1)+Code2)
        F.NT(1) = NT1;
        F.NT(2) = NT2;
        F.Filename = E.Filename;
        zDisplayNT(F,1:2,VP);
        view(2)
      end
    end
  end

pause

end


return

load PairExemplars

[s,t] = size(Exemplars);

%for r = 1:s,                              % row of the table
%  for pc = 1:t,                           % paircode of first pair
%    if ~isempty(Exemplar(r,pc).Filename),
%      N1 = Exemplar(r,pc).NT1;
%      N2 = Exemplar(r,pc).NT2;
%      for q = 1:s,                        % row of table
%        for qc = 1:t,                     % paircode of second pair
%          if ~isempty(Exemplar(q,qc).Filename),
%            M1 = Exemplar(r,pc).NT1;
%            M2 = Exemplar(r,pc).NT2;
%            d = zIsoDiscrepancy(N1,N2,M1,M2);
