
Filenames = 'Nonredundant_2008_06_06_list';

% ----------------------------------------- Read PDB structure names and lists

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end

% FullList = FullList(1:50);

% --------------------------------- Look up information on these structures

load PDBInfo

i = [];

for f = 1:length(FullList),
  pp = find(ismember(upper(t(:,1)),upper(FullList{f})));
  if ~isempty(pp),
    i = [i pp(1)];   
  else
    fprintf('No information on file %s\n', FullList{f});
  end
end
t = t(i,:);                       % omit structures with no information
n = n(i,:);           

i = find(n(:,3) > 0);              % omit structures with no basepairs
t = t(i,:);                       % omit structures with no information
n = n(i,:);           

F = length(i);                    % number of files

[y,i] = sort(n(:,2));             % sort by number of nucleotides

t = t(i,:);
n = n(i,:);

% ---------------------------------- List information on these files

for f = 1:F,
  g{f,1} = t{f,1};                % filename
  g{f,2} = n(f,1);                % resolution
  g{f,3} = n(f,2);                % number of nucleotides
  g{f,4} = n(f,3);                % number of basepairs
  g{f,5} = t{f,2};                % descriptor
  g{f,6} = t{f,3};                % experimental method
  g{f,7} = t{f,4};                % release date
  g{f,8} = t{f,5};                % author
  g{f,9} = t{f,6};                % keywords
  g{f,10} = t{f,8};               % biological source
end

g = [{'Filename','Resolution','NumNucleotides','NumBasepairs','Descriptor','ExpMethod','ReleaseDate','Author','Keywords','Source'}; g];

xlswrite('ReportOnList.xls',g)

% ----------------------------------- Histogram of number of nucleotides

Edges = [10*(0:20) 100*(3:10) 500*(3:6)];
N = histc(n(:,2),Edges);
figure(1)
clf
%bar(N,'histc')
bar(0:(length(Edges)-1),N,0.8,'histc')
for i = 1:length(Edges),
  Lab{i} = num2str(Edges(i));
end
set(gca,'XTick',(0:(length(Edges)-1)))
set(gca,'XTickLabel',Lab,'FontSize',8)
xlabel('Number of nucleotides. Note the irregular scale.');
ylabel('Number of distinct instances within each range');
title('Illustration of the number of nucleotides in the distinct instances');

% ----------------------------------- Histogram of number of basepairs

Edges = [0 1 10*(1:20) 100*(3:15)];
N = histc(n(:,3),Edges);
figure(2)
clf
%bar(N,'histc')
bar(0:(length(Edges)-1),N,0.8,'histc')
for i = 1:length(Edges),
  Lab{i} = num2str(Edges(i));
end
set(gca,'XTick',[0.5 (1:(length(Edges)-1))])
set(gca,'XTickLabel',Lab,'FontSize',8)
xlabel('Number of basepairs. Note the irregular scale.');
ylabel('Number of distinct instances within each range');
title('Illustration of the number of basepairs in the distinct instances');

% ----------------------------------- Scatterplot basepairs vs nucleotides

% Note: several structures have the same number of basepairs and nucleotides,
% but they appear as a single point on the graph.  We could modify the size
% of the symbol according to how many instances there are.  That will work 
% better than shifting them randomly.


figure(3)
clf
loglog(n(:,2),n(:,3),'.')
xlabel('Number of nucleotides, minimum 1');
ylabel('Number of basepairs, minimum 1');
title('Non-redundant dataset summary');

FN = '1FFK';
Text = '23S rRNA ';

i = find(ismember(t(:,1),FN));
text(n(i,2),n(i,3),Text,'HorizontalAlignment','Right');

i = find(ismember(t(:,1),'2AVY'));
text(n(i,2),n(i,3),'16S rRNA ','HorizontalAlignment','Right');

i = find(ismember(t(:,1),'1X8W'));
text(n(i,2),n(i,3),'Tetrahymena Ribozyme ','HorizontalAlignment','Right');

i = find(ismember(t(:,1),'3BWP'));
text(n(i,2),n(i,3),'Group II intron  ','HorizontalAlignment','Right');

i = find(ismember(t(:,1),'1YFG'));
text(n(i,2),n(i,3),'Yeast initiator tRNA       ','HorizontalAlignment','Right');

