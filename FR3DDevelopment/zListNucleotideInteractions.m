
function [void] = zListNucleotideInteractions(File,NTList)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1 || isempty(NTList),
  NTList = 1:File.NumNT;                  % display them all
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

for a = 1:length(Indices),
  i = Indices(a);
  e = File.Edge(i,:);                         % interactions it makes
  j = find((abs(e) > 0) .* (abs(e) < 100));   % NTs i interacts with
  [y,k] = sort(abs(e(j)));                    % cWW first
  j = j(k);

  fprintf('%4s %1s%4s interactions: ', File.Filename, File.NT(i).Base, File.NT(i).Number);
  for b = 1:length(j),
    if b > 1,
      fprintf('|');
    else
      fprintf(' ');
    end
    fprintf('%5s %1s%4s ', zEdgeText(e(j(b)),1,File.NT(i).Code,File.NT(j(b)).Code), File.NT(j(b)).Base, File.NT(j(b)).Number);
  end

  bp = File.BasePhosphate(i,:);

  j = find((abs(bp) > 0) .* (abs(bp) < 100));   % NTs i interacts with
  [y,k] = sort(abs(bp(j)));                    % put in standard order
  j = j(k);

  for b = 1:length(j),
    t = zBasePhosphateText(bp(j(b)));
    if t(2) ~= '0',
      fprintf('#');
      fprintf('%5s %1s%4s ', t, File.NT(j(b)).Base, File.NT(j(b)).Number);
    end
  end

  bp = File.BasePhosphate(:,i);

  j = find((abs(bp) > 0) .* (abs(bp) < 100));   % NTs i interacts with
  [y,k] = sort(abs(bp(j)));                    % put in standard order
  j = j(k);

  for b = 1:length(j),
    t = zBasePhosphateText(-bp(j(b)));
    if t(2) ~= '0',
      fprintf('#');
      fprintf('%5s %1s%4s ', t, File.NT(j(b)).Base, File.NT(j(b)).Number);
    end
  end

  fprintf('\n');
end

    