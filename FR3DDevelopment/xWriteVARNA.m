
function [T] = xWriteVARNA(F,NTList,Verbose)

if nargin < 3,
  Verbose = 0;
end

% ---------------------------------------------- Determine type of file

if isfield(F,'OKCodes'),
  Type = 'Q';                                  % query
else
  Type = 'F';                                  % file
end

% ----------------------------------------------- Interpret NTList, F

if Type == 'Q',                                % query
  if isempty(NTList),
    Indices = 1:F.NumNT;
  else
    Indices = NTList;
  end
else
  if strcmp(class(F),'char'),
    Filename = F;
    F = zGetNTData(Filename,0);
  end

  if nargin == 1 || isempty(NTList),
    NTList = 1:length(File.NT);                  % display them all
  end

  % if NTList is a cell array of numbers, look up the indices

  if strcmp(class(NTList),'char'),
    NTList = {NTList};
  end

  if strcmp(class(NTList),'cell'),
    Indices = zIndexLookup(F,NTList);
  else
    Indices = NTList;
  end
end

N = length(Indices);

if Type == 'Q',
  for i = 1:N,
    for j = 1:N,
      if ~isempty(F.EdgeNums{i,j}),
        [y,a] = min(abs(F.EdgeNums{i,j}));       % location of best option
        if y < 13,                               % a basepair
          BP(i,j) = F.EdgeNums{i,j}(a);
        else
          BP(i,j) = 0;
        end
      else
        BP(i,j) = 0;
      end
    end
  end
else
  BP = F.Edge(Indices,Indices);                  % extract relevant interactions
  BP = BP .* (abs(BP) > 0) .* (abs(BP) < 13);
end

% ---------------------------------------------- Codes for sequence signatures

Code(8) = 'A';
Code(4) = 'C';
Code(2) = 'G';
Code(1) = 'U';
Code(12) = 'M';
Code(10) = 'R';
Code(9) = 'W';
Code(6) = 'S';
Code(5) = 'Y';
Code(3) = 'K';
Code(14) = 'V';
Code(13) = 'H';
Code(11) = 'D';
Code(7) = 'B';
Code(15) = 'N';

% ----------------------------------------------- Produce text

c = 1;
T{c} = '<html>';
c = c + 1;
T{c} = '<applet code="VARNA.class"';
c = c + 1;
T{c} = '	codebase="bin"';
c = c + 1;
T{c} = '	archive="VARNA.jar"';
c = c + 1;
T{c} = '	width=300" height="300">';
c = c + 1;
T{c} = '	<param name="sequenceDBN"  value="';

if Type == 'Q',
  for i = 1:N,
    j = F.OKCodes{i} * [8 4 2 1]';
    T{c} = [T{c} Code(j)];
  end
else
  for i = 1:N,
    T{c} = [T{c} F.NT(Indices(i)).Base];
  end
end

T{c} = [T{c} '" />'];
c = c + 1;

T{c} = '	<param name="structureDBN" value="';
for i = 1:N,
  T{c} = [T{c} '.'];
end
T{c} = [T{c} '" />'];
c = c + 1;

if sum(sum(abs(BP))) > 0,
  T{c} = '	<param name="auxBPs" value="';
  c = c + 1;
  for i = 1:N,
    for j = (i+1):N,
      if abs(BP(i,j)) > 0,                     % basepair here
        t = lower(zEdgeText(BP(i,j)));         % type of pair
        if t(1) == 'c',
          g = 'cis';
        else
          g = 'trans';
        end
        if t(2) == 'w',
          e = 'wc';
        else
          e = t(2);
        end
        if t(3) == 'w',
          f = 'wc';
        else
          f = t(3);
        end
        T{c} = sprintf('(%d,%d):edge5=%s,edge3=%s,stericity=%s;',i,j,e,f,g);
        c = c + 1;
      end
    end
  end
end

T{c} = '	<param name="rotation" value="-40" />';
c = c + 1;
T{c} = '	<param name="bpStyle" value="lw" />';
c = c + 1;
T{c} = '	<param name="background" value="#FFFFFF" />';
c = c + 1;
T{c} = '</applet>';
c = c + 1;
T{c} = '</html>';

if Verbose > 0,
  for c = 1:length(T),
    fprintf('%s\n', T{c});
  end
end
