% zPhosphateLocations reads a secondary structure annotation from an Excel spreadsheet, then tabulates the number of local and long-range BPh interactions according to the secondary structure element in which the base and the phosphate are found


if 0 > 1,
 Filenames = {'2AW4','2J01','2AVY','1J5E'};
 Filenames = {'2J01','1J5E'};                % T. Th. 16S, 23S
 Filenames = {'2AVY','2AW4'};                % E. coli 16S, 23S

 File = zAddNTData(Filenames);

 for f = 1:length(File),
  %--------  make the next line specific to the structure file of interest
  [n,t,r] = xlsread(['Annotation' filesep File(f).Filename '_Annotation.xls']);

  % -------- Annotate the structure file with secondary structure information
  % -------- This technique is slow.  Faster to go line by line, checking Num
  for a = 1:length(n(:,1)),
    if strcmp(class(r{a,2}),'double'),
      num = num2str(r{a,2});
    else
      num = r{a,2};
    end
    i = zIndexLookup(File(f),num,r{a,3});
    File(f).NT(i(1)).Hierarchy{1} = '';
    File(f).NT(i(1)).Hierarchy{2} = '';
    File(f).NT(i(1)).Hierarchy{3} = upper(r{a,4});
    File(f).NT(i(1)).Hierarchy{4} = '';
  end

  % ------- List nucleotides that might be in the wrong secondary str element

  fprintf('Listing for %4s\n', File(f).Filename);
  for i = 1:length(File(f).NT),
   if ~isempty(File(f).NT(i).Hierarchy),
    if isempty(strfind(File(f).NT(i).Hierarchy{3},'HL')) && ...
       ~isempty(strfind(File(f).NT(i).Hierarchy{3},'H')) && ...
       length(find(fix(abs(File(f).Edge(i,:))) == 1)) == 0,
      fprintf('File %4s nucleotide %s%4s is in %4s but makes no cWW\n', ...
        File(f).Filename, File(f).NT(i).Base, File(f).NT(i).Number, ...
        File(f).NT(i).Hierarchy{3});
    end
   end
  end
 end
end


c = 1;
clear L                           % place to accumulate interactions
for f = 1:length(File),

  % -------- Find the base-phosphate interactions

  [i,j] = find((File(f).BasePhosphate > 0) .* (File(f).BasePhosphate < 20));

  k = find(i~=j);                             % keep non-self interactions
  i = i(k);
  j = j(k);

  for k = 1:length(i),
    if length(File(f).NT(i(k)).Hierarchy) > 0,
      L{c,1} = [File(f).NT(i(k)).Hierarchy{3} num2str(f)];     % element containing the base
    else
      L{c,1} = ''; 
    end
    if length(File(f).NT(j(k)).Hierarchy) > 0,
      L{c,2} = [File(f).NT(j(k)).Hierarchy{3} num2str(f)];     % element containing the phos
    else
      L{c,2} = '';
    end
    L{c,3} = min(File(f).Range(i(k),j(k)),abs(i(k)-j(k))-1); % interaction range
    L{c,4} = File(f).NT(i(k)).Base;
    L{c,5} = File(f).NT(i(k)).Number;
    L{c,6} = File(f).NT(j(k)).Base;
    L{c,7} = File(f).NT(j(k)).Number;
    L{c,8} = zBasePhosphateText(File(f).BasePhosphate(i(k),j(k)));
    L{c,9} = zEdgeText(File(f).Edge(i(k),j(k)));
    L{c,10} = f;
    L{c,11} = File(f).Crossing(i(k),j(k));   % number of cWW's crossed
    L{c,12} = double(strcmp(L{c,1},L{c,2}));
    c = c + 1;
  end
end

[y,i] = sort(cat(1,L{:,11}));
[y,i] = sort(L(:,1));
L = L(i,:);

%L(:,[4 5 1 6 7 2 8 3 9 11 10])

fprintf('Base-phosphate interactions between secondary structure elements.\n');
fprintf('Files analyzed are ');
for f = 1:length(File),
  fprintf('%s ', File(f).Filename);
end
fprintf('\n');

for y = 1:2,                                        % two tables, short/long
  switch y
  case 1, 
    i = find(cat(1,L{:,11}) <  2);                  % short range interactions
    ii = find(cat(1,L{:,12}) > 0);                  % same element name
    i = union(i,ii);
    fprintf('Short range phosphate interactions\n');
  case 2, 
    i = find((cat(1,L{:,11}) >= 2).*(cat(1,L{:,12}) == 0)); % long range

    fprintf('Long range phosphate interactions\n');
  end

  M = L(i,:);                                       % focus on these

  Count = zeros(10,10);                             % accumulate counts

  clear MM
  for j = 1:length(i),
    a = 0;
    b = 0;
    if ~isempty(M{j,1}) && ~isempty(M{j,2}),
      a = zAnnotationType(upper(M{j,1}));
      b = zAnnotationType(upper(M{j,2}));
      MM(j,1) = a;
      MM(j,2) = b;
      Count(a,b) = Count(a,b) + 1;                  % accumulate counts

if y==2 && a==5 && b==5,
%  M(j,:)
end

if y==1 && a==2 && b==2,
  M(j,:)
end

    end
  end

  [z,k] = sortrows(MM,[1 2]);

  Type = {'Hairpin','Helix','Internal','Bulge','Junction','End'};

  h = [1 2 3 5];

  fprintf('     Phosphate -->  ');
  for b = h,
    fprintf('%10s ', Type{b});
  end
  fprintf('     Total\n');  

  for a = h,
    fprintf('         %10s ', Type{a});
    for b = h,
      fprintf('%10d ', Count(a,b));
    end
    fprintf('%10d\n', sum(Count(a,:)));
  end

  fprintf('              Total ');
  for b = h,
    fprintf('%10d ', sum(Count(:,b)));
  end
  fprintf('%10d\n', sum(sum(Count)));

end

% ------------------------------- How many hairpins contain a BPh? etc.


% HL  -> 1
% H   -> 2
% IL  -> 3
% BL  -> 3   % treat it the same as an internal loop
% J   -> 5
% 3_p -> 6
% 5_p -> 7

J = unique(L(:,1));                       % list each element once
C = zeros(1,5);

for j = 1:length(J),
  if ~isempty(strfind(J{j},'HL')),
    C(1) = C(1) + 1;
  elseif ~isempty(strfind(J{j},'H')),
    C(2) = C(2) + 1;
  elseif ~isempty(strfind(J{j},'IL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'BL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'J')),
    C(5) = C(5) + 1;
  elseif ~isempty(strfind(J{j},'3_p')),
    C(6) = C(6) + 1;
  elseif ~isempty(strfind(J{j},'5_p')),
    C(7) = C(7) + 1;
  end
end

D = C;

OK = zeros(length(L(:,1)));

for i = 1:length(L(:,1)),
  if strcmp(L{i,1},L{i,2}) && ~isempty(L{i,1}),    % internal to an element
    OK(i) = 1;
  end
end

J = L(find(OK),1);

C = zeros(1,5);

for j = 1:length(J),
  if ~isempty(strfind(J{j},'HL')),
    C(1) = C(1) + 1;
  elseif ~isempty(strfind(J{j},'H')),
J{j}
    C(2) = C(2) + 1;
  elseif ~isempty(strfind(J{j},'IL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'BL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'J')),
    C(5) = C(5) + 1;
  elseif ~isempty(strfind(J{j},'3_p')),
    C(6) = C(6) + 1;
  elseif ~isempty(strfind(J{j},'5_p')),
    C(7) = C(7) + 1;
  end
end

E = C;

J = unique(L(find(OK),1));

C = zeros(1,5);

for j = 1:length(J),
  if ~isempty(strfind(J{j},'HL')),
    C(1) = C(1) + 1;
  elseif ~isempty(strfind(J{j},'H')),
    C(2) = C(2) + 1;
  elseif ~isempty(strfind(J{j},'IL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'BL')),
    C(3) = C(3) + 1;
  elseif ~isempty(strfind(J{j},'J')),
    C(5) = C(5) + 1;
  elseif ~isempty(strfind(J{j},'3_p')),
    C(6) = C(6) + 1;
  elseif ~isempty(strfind(J{j},'5_p')),
    C(7) = C(7) + 1;
  end
end

fprintf('Number of secondary structure elements containing at least one BPh interaction\n');

fprintf('Hairpins     %3d      %3d      %3d\n', C(1), D(1), E(1));
fprintf('Helices      %3d      %3d      %3d\n', C(2), D(2), E(2));
fprintf('Internal     %3d      %3d      %3d\n', C(3), D(3), E(3));
fprintf('Junctions    %3d      %3d      %3d\n', C(5), D(5), E(5));
