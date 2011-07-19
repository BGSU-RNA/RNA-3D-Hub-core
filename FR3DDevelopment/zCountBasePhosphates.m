% zCountBasePhosphates(File) tabulates the base-phosphate interactions in File
% It formats the counts by geometric family and, within each family, by basepair.
% It writes an Excel sheet with a standard format.
% CountFilename could be 'Basepair_counts_NonRedundant_2008_02_21_list.xls'
% If omitted, the filenames in File are used to construct the filename.

% File = zAddNTData('NonRedundant_2008_02_21_list');
% xAnnotate
% F    = zResolutionFilter(File,2.5);
% zCountBasePhosphates(F,'BPh_counts');

function [Count,CCount] = zCountBasePhosphates(File,CountFilename)

if nargin == 1,
  CountFilename = 'BasePhosphate_counts';
  for f = 1:length(File),
    CountFilename = [CountFilename '_' File(f).Filename];
  end
  CountFilename = [CountFilename '.xls'];
end

MaxCat = 25;                                  % maximum category number

CI  = [];                                     % count interactions
Data = [];

for f = 1:length(File),
  B = File(f).BasePhosphate;                  % BPh interactions
  for i = 1:length(B(:,1)),
    B(i,i) = 0;                               % eliminate self interactions
  end
  B = B .* (B > 0) .* (B < MaxCat);           % only certain interactions
  [i,j,k] = find(B);
  Code1 = cat(1,File(f).NT(i).Code);
  Code2 = cat(1,File(f).NT(j).Code);
  
  CI  = [CI; [k Code1 Code2]];

  r = File(f).Info.Resolution;
  if isempty(r),
    r = NaN;
  end

  E = fix(abs(File(f).Edge));  
  ncww = sum(sum( (E > 1) .* (E < 13)));       % number of non-cWW
  cww = sum(sum((E == 1)));       % number of cWW

  NewData = full([f sum(sum(B > 0)) length(File(f).NT) r cat(2,File(f).Motifs(1:6).Count) ncww/length(File(f).NT) cww/length(File(f).NT)]);
  Data = [Data; NewData];

end

Lab{1} = 'Filename';
Lab{2} = 'Number of BPh';
Lab{3} = 'Number of nucleotides';
Lab{4} = 'Resolution';
Lab{5} = 'Number of pseudoknots';
Lab{6} = 'Number of 6-nucleotide helices';
Lab{7} = 'Number of hairpins';
Lab{8} = 'Number of internal loops';
Lab{9} = 'Number of bulge loops';
Lab{10} = 'Number of junction strands';
Lab{11} = 'non-cWW per nucleotide';
Lab{12} = 'cWW per nucleotide';

% ---------------------------------- plot num BPh against molecule size

Data(:,10) = Data(:,10) - 2*Data(:,8) - Data(:,7);

Data = [Data Data(:,2)./Data(:,3)];
[s,t] = size(Data);

figure(1)
clf
i = find(Data(:,4) > 2.5);
semilogx(Data(i,3),Data(i,t),'r.')
hold on
i = find(Data(:,4) <= 2.5);
semilogx(Data(i,3),Data(i,t),'b.')
xlabel('Number of nucleotides');
ylabel('BPh interaction per nucleotide');

figure(2)
clf
i = find(Data(:,4) > 2.5);
plot(Data(i,7),Data(i,t),'r.');
hold on
i = find(Data(:,4) <= 2.5);
plot(Data(i,7),Data(i,t),'b.');
xlabel('Number of hairpins');
ylabel('BPh interactions per nucleotide');

figure(3)
clf
i = find(Data(:,4) > 2.5);
plot(Data(i,4),Data(i,t),'r.');
hold on
i = find(Data(:,4) <= 2.5);
plot(Data(i,4),Data(i,t),'b.');
xlabel('Resolution');
ylabel('BPh interactions per nucleotide');

%semilogx(Data(:,3),Data(:,2),'.')
%plot(Data(:,3),Data(:,2),'.')
%loglog(Data(:,3),Data(:,2),'.')

figure(4)
clf
i = find(Data(:,4) > 2.5);
plot(Data(i,t-1)/2,Data(i,t),'r.');
hold on
i = find(Data(:,4) <= 2.5);
plot(Data(i,t-1)/2,Data(i,t),'b.');
xlabel('cWW per nucleotide');
ylabel('BPh interactions per nucleotide');
saveas(gcf,'Bph per NT versus cWW per NT.png','png');

figure(5)
clf
i = find(Data(:,4) > 2.5);
plot(Data(i,7)./Data(i,3),Data(i,t),'r.');
hold on
i = find(Data(:,4) <= 2.5);
plot(Data(i,7)./Data(i,3),Data(i,t),'b.');
xlabel('Hairpins per nucleotide');
ylabel('BPh interactions per nucleotide');

figure(6)
clf
i = find(Data(:,4) > 2.5);
plot(Data(i,t-2)/2,Data(i,t),'r.');
hold on
i = find(Data(:,4) <= 2.5);
plot(Data(i,t-2)/2,Data(i,t),'b.');
xlabel('non-cWW per nucleotide');
ylabel('BPh interactions per nucleotide');
saveas(gcf,'Bph per NT versus cWW per NT.png','png');

i = find(Data(:,t-1) > 0);            % files with cWW's

Vars = [5:10];
Vars = [9 10 11 12];
Vars = [11 12];
Vars = [5:10 12];
Vars = [11 12];

[B,BINT,R,RINT,STATS] = regress(Data(i,t),[Data(i,Vars) ones(length(i),1)]);

STATS
for v = 1:length(Vars),
  vv = Vars(v);
  fprintf('Variable %2d is %30s interval %7.4f to %7.4f\n', vv, Lab{vv}, BINT(v,1), BINT(v,2));
end
v = length(Vars) + 1;
fprintf('Constant term           interval %7.4f to %7.4f\n', BINT(v,1), BINT(v,2));

[B,BINT,R,RINT,STATS] = regress(Data(i,t),[Data(i,11)./Data(i,12) ones(length(i),1)]);


figure(7)
clf
scatter3(Data(:,11),Data(:,12),Data(:,13),10,(Data(:,3)>50),'filled');
%scatter3(Data(:,11),Data(:,12),Data(:,13),10,log(Data(:,3)),'filled');
xlabel('Non-cWW per NT');
ylabel('cWW per NT');
zlabel('BPh frequency');
axis([0 1 0 1 0 0.25]);
rotate3d on
caxis([-0.2 1.2])


Data = sortrows(Data,[-t -3]);

for i = 1:length(Data(:,1)),
  f = Data(i,1);
  fprintf('File %4s NTs %4d BPh %4d BPhPercentage %7.2f cWW/NT %7.2f ', File(f).Filename, Data(i,3), Data(i,2), Data(i,t), Data(i,t-1)/2);

  fprintf('Resolution: %4.1f ', Data(i,4));

  fprintf('%3d ', Data(i,5:10));

  fprintf('%s\n', File(f).Info.Descriptor);
end

Count = zeros(4,4,MaxCat);

for i = 1:length(CI(:,1)),
  Count(CI(i,2),CI(i,3),CI(i,1)) = Count(CI(i,2),CI(i,3),CI(i,1)) + 1;
end

BPCat = [2 6 7 0 6 7 8 9 0 1 3 4 5 0 5 9 0 7 4];  % updated 8-19-2008

% ---------------------------------------- start producing output

row = 1;
clear T
Letters = 'ACGU';

% ---------------------------------------- counts by interacting base

NC = [];

for i = 1:4,                               % base making the base-phosphate
  for j = 0:9,                             % interaction being made
    NC(i,j+1) = sum(sum(Count(i,:,find(BPCat == j))));
  end
end

T{row,1} = 'Base nucleotide';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
  T{row+i,1} = Letters(i);
  for j = 0:9,
    T{row+i,2+j} = NC(i,j+1);
  end
  T{row+i,12} = sum(NC(i,2:10));
  T{row+i,13} = sum(NC(i,:));
end

% ---------------------------------------- counts by phosphate identity

NC = [];

row = row + 6;

for i = 1:4,                               % base making the base-phosphate
  for j = 0:9,                             % interaction being made
    NC(i,j+1) = sum(sum(Count(:,i,find(BPCat == j))));
  end
end

T{row,1} = 'Phosphate nucleotide';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
  T{row+i,1} = Letters(i);
  for j = 0:9,
    T{row+i,2+j} = NC(i,j+1);
  end
  T{row+i,12} = sum(NC(i,2:10));
  T{row+i,13} = sum(NC(i,:));
end

% ---------------------------------------- counts by pair of bases

row = row + 6;
NC = [];

for i = 1:4,                               % base making the base-phosphate
 for k = 1:4,
  for j = 0:9,                             % interaction being made
    NC(i,k,j+1) = sum(Count(i,k,find(BPCat == j)));
  end
 end
end

T{row,1} = 'Pair of bases';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
 for k = 1:4,
  p = 4*(i-1) + k;
  T{row+p,1} = Letters([i k]);
  for j = 0:9,
    T{row+p,2+j} = NC(i,k,j+1);
  end
  T{row+p,12} = sum(NC(i,k,2:10));
  T{row+p,13} = sum(NC(i,k,:));
 end
end

% --------------------------------------- write to spreadsheet

if ~isempty(CountFilename),
  xlswrite(CountFilename,T);
end


