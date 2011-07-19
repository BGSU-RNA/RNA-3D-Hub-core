% zConservationPrediction

% Read Jesse's dataset on E.coli with number of interactions and potential base-protein interactions

[n,t] = xlsread('Bases_from_Ec_and_Interactions_10_21-1.xls');

t = t(3:end,:);                               % remove first two rows

Filenames = {'2AVY','2AW4'};

 File = zAddNTData(Filenames);
% File = zAttachAlignment(File,1);
% xAnnotate

Data = [];

for f = 1:length(File),
  BP(f).BP = File(f).BasePhosphate;
  for nn = 1:File(f).NumNT,
    BP(f).BP(nn,nn) = 0;                      % no self interactions
  end

  E = abs(fix(File(f).Edge));
  E = (E > 0) .* (E < 25);               % basepairing indicator

  File(f).BaseInter = E + E^2;              % direct or indirect connection
end

for r = 1:length(n(:,1)),
  f = find(ismember(lower(Filenames),lower(t{r,1})));
  % fprintf('%d %s  %s\n', f, Filenames{f}', t{r,1})
  Data(r,1) = f;
  switch f,
    case 1, i = zIndexLookup(File(f),num2str(n(r,1)),'A');
    case 2, i = zIndexLookup(File(f),num2str(n(r,1)),'B');
  end
  i = i(1);
  Data(r,2) = i;                                % index of this nucleotide
  Data(r,17) = File(f).NT(i).Code;              % which base this is

  Data(r,8) = sum((abs(File(f).Edge(i,:)) > 20) .* (abs(File(f).Edge(i,:)) < 24));

  Data(r,10) = sum(((BP(f).BP(i,:)) > 0) .* ((BP(f).BP(i,:)) < 30));
  Data(r,12) = sum(((BP(f).BP(:,i)) > 0) .* ((BP(f).BP(:,i)) < 30));

  Data(r,15) = 0;
  Data(r,16) = 0;

  if Data(r,12) > 0,
    j = find(((BP(f).BP(:,i)) > 0) .* ((BP(f).BP(:,i)) < 30));
    for k = 1:length(j),
      e = abs(fix(File(f).BaseInter(i,j(k))));         % edge between i, j(k)
      b = abs(fix(BP(f).BP(i,j(k))));             % BPh with i and j(k)
      if (e > 0),                                 % i also pairs with j(k)
        Data(r,16) = Data(r,16) + 1;
      else
        Data(r,15) = Data(r,15) + 1;
      end
    end
  end

  Data(r,13) = sum((fix(abs(File(f).Edge(i,:))) ==1));
  Data(r,14) = sum((abs(File(f).Edge(i,:)) >= 2) .* (abs(File(f).Edge(i,:)) < 15));

  counts = n(r,4:7);
  percs  = 0.00001 + counts / sum(counts);        % number of A, C, G, U

  Data(r,11) = - percs * log2(percs)';            % entropy of the distribution
  Data(r,3) = 100*counts(File(f).NT(i).Code)/sum(counts); % conservation

  File(f).Conservation(i) = Data(r,3);

  if isnan(Data(r,11)),
    counts
    Data(r,3)
  end
end

cp = fix(mean(Data(:,3)));
i = find(Data(:,8) == 0);                      % NT with no stacking partners
cp = fix(mean(Data(i,3)));

for r = 1:length(n(:,1)),
  f = Data(r,1);
  i = Data(r,2);
  e = abs(fix(File(f).Edge(i,:)));                  % interaction list
  j = find((e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    Data(r,18) = mean(c)/100;             % conservation of stacking partners
  else
    Data(r,18) = cp/100;
    if i == 1,
      Data(r,18) = File(f).Conservation(i+1)/100;
    elseif i == File(f).NumNT,
      Data(r,18) = File(f).Conservation(i-1)/100;
    else
      Data(r,18) = mean(File(f).Conservation([i-1 i+1]))/100;
    end
  end
end

sum(Data(:,15:end))                   % see how many phosphate acceptors 

% Data(:,1)                           % file number
% Data(:,2)                           % index of nucleotide in the file
% Data(:,3)                           % my calculation of conservation %

Data(:,4) = n(:,12);                  % Jesse's count of #cWW
Data(:,5) = n(:,13);                  % Jesse's count of #non-cWW
Data(:,6) = n(:,14);                  % Jesse's count of #BPh
Data(:,7) = n(:,15);                  % Potential base-protein

% Data(:,8)                           % my count of stacking interactions

Data(:,9) = ones(size(Data(:,8)));    % constant term

% Data(:,10)                          % my count of number of times donor
% Data(:,11)                          % entropy of the distribution
% Data(:,12)                          % my count of number of times acceptor
% Data(:,13)                          % my count of #cWW
% Data(:,14)                          % my count of #non-cWW
% Data(:,15)                          % acceptor, no connection to donor
% Data(:,16)                          % acceptor, connection to donor
% Data(:,17)                          % base code
% Data(:,18)                          % conservation of stacking partners
% Data(:,19)                          % predicted conservation percentage

for r = 1:length(Data(:,1)),
  for c = 1:length(Data(1,:)),
    if isnan(Data(r,c)),
      Data(r,c) = 0;
    end
  end
end

% Compare Jesse's listing to mine for cWW and non-cWW
%[Data(:,4) Data(:,13) (Data(:,4)-Data(:,13))]
%[Data(:,5) Data(:,14) (Data(:,5)-Data(:,14))]
%[Data(:,6) Data(:,10) (Data(:,6)-Data(:,10))]

Response = Data(:,11);                % entropy of the base distribution
Range    = 0.05+[0:0.1:max(Response)];
Resp     = 'Entropy of base position in sequences';

Response = Data(:,3);                 % conservation percentage
Range    = 1+[0:2:98];
Resp     = 'Conservation percentage of base in sequences';


VarName{1} = 'file number';
VarName{2} = 'index of nucleotide';
VarName{3} = 'my calculation of conservation percentage';
VarName{4} = 'number of cWW pairs';
VarName{5} = 'number of non-cWW pairs';
VarName{6} = 'number of base-phosphate interactions';
VarName{7} = 'near a protein (1) or not (0)';
VarName{8} = 'number of stacking interactions';
VarName{9} = '(constant term)';
VarName{10} = 'number of BPh interactions in which it is the H-bond donor';
VarName{11} = 'entropy of base distribution';
VarName{12} = 'number of BPh interactions in which it is the phosphate acceptor';
VarName{13} = 'FR3D cWW pair';
VarName{14} = 'FR3D non-cWW pair';
VarName{15} = 'number of times phosphate acceptor but no chain to base donor';
VarName{16} = 'number of times phosphate acceptor with chain to base donor';
VarName{17} = 'base code';
VarName{18} = 'conservation percentage of stacking partners or neighboring bases, divided by 100';

Variables = [4 5 10 7 8 12 9];  % put the constant term at the end!
Variables = [4 5 10 7 8 15 16 9];  % put the constant term at the end!
Variables = [4 5 10 7 8 15 9];  % put the constant term at the end!
Variables = [4 5 10 7 8 16 9];  % put the constant term at the end!
Variables = [18 5 10 16 7 15 4 9];  % put the constant term at the end!
Variables = [18 5 10 16 7 15 4 9];  % put the constant term at the end!
Variables = [18 5 10 7 12 4 9];  % put the constant term at the end!

StData = Data(:,Variables);

max(StData);

[B,BINT,R,RINT,STATS] = regress(Response,StData);

Data(:,19) = StData * B;

fprintf('Average %s is %7.4f\n', Resp, mean(Response));

fprintf('Linear regression model for all nucleotides:\n');
fprintf('%s = \n', Resp);
fprintf('%7.1f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.1f * %15s\n', abs(B(v)), VarName{Variables(v)});
end


BINT

Variables = [5 10 7 12 4 18 9];  % put the constant term at the end!
for v = 1:(length(Variables)-1),
  V = [Variables(1:v) 9];
  StData = Data(:,V);
  max(StData);
  [B,BINT,R,RINT,STATS] = regress(Response,StData);
  V
  STATS
end



Variables = Variables(2:end);
StData = StData(:,2:end);


zConservationPredictionExcel


break

if 0 > 1,
r = find(Data(:,6) ==3);                   % NT with three BPh
for j = 1:length(r),
  f = Data(r(j),1);
  n = Data(r(j),2);
  fprintf('%d %d %s\n', r(j), n, File(f).NT(n).Number);
  k = find(File(f).BasePhosphate(n,:));
  zShowInteractionTable(File(f),[n k]);
  VP.Sugar = 1;
  clf
  zDisplayNT(File(f),[n k],VP);
  pause
end

r = find(Data(:,8) == 4);                   % NT with four stacks
for j = 1:length(r),
  f = Data(r(j),1);
  n = Data(r(j),2);
  fprintf('%d %d %s\n', r(j), n, File(f).NT(n).Number);
  k = find(File(f).Edge(n,:));
  zShowInteractionTable(File(f),[n k]);
  VP.Sugar = 1;
  clf
  zDisplayNT(File(f),[n k],VP);
  pause
end

r = find(Data(:,5) == 4);                   % NTs with 4 non-cWW
r = find(Data(:,5) == 3);                   % NTs with 3 non-cWW
for j = 1:length(r),
  f = Data(r(j),1);
  n = Data(r(j),2);
  fprintf('%d %d %s\n', r(j), n, File(f).NT(n).Number);
  k = find(File(f).Edge(n,:));
  zShowInteractionTable(File(f),[n k]);
  VP.Sugar = 1;
  clf
  zDisplayNT(File(f),[n k],VP);
  pause
end
end

i = find(Data(:,4) > 0);             % cWW pairs
[B,BINT,R,RINT,STATS] = regress(Response(i),StData(i,:));
fprintf('Linear regression model for nucleotides making a cWW pair:\n');
fprintf('%s = \n', Resp);
fprintf('%7.4f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.4f * %15s\n', abs(B(v)), VarName{Variables(v)});
end

i = find(Data(:,4) == 0);             % not making a cWW pair
[B,BINT,R,RINT,STATS] = regress(Response(i),StData(i,:));
fprintf('Linear regression model for nucleotides not making a cWW pair:\n');
fprintf('%s = \n', Resp);
fprintf('%7.4f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.4f * %15s\n', abs(B(v)), VarName{Variables(v)});
end

% break

fn = 1;

figure(fn)
clf
subplot(1,2,1)
hist(Response,Range)
Text = 'Histogram of nucleotide conservation - all nucleotides'
title(Text);
subplot(1,2,2)
boxplot(Response);
set(gcf,'Renderer','painters');    % makes nice PDF files
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
c = 4;
i = find(Data(:,c) > 0);
subplot(2,2,1)
hist(Response(i),Range)
Text = 'Nucleotide conservation - nucleotides making cWW';
title(Text);
subplot(2,2,2)
boxplot(Response(i));
subplot(2,2,3)
i = find(Data(:,c) == 0);
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides not making cWW');
subplot(2,2,4)
boxplot(Response(i));
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
c = 5;
i = find(Data(:,c) > 0);
subplot(2,2,1)
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides making non-cWW');
subplot(2,2,2)
boxplot(Response(i));
subplot(2,2,3)
i = find(Data(:,c) == 0);
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides not making non-cWW');
subplot(2,2,4)
boxplot(Response(i));
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
c = 8;
i = find(Data(:,c) > 0);
subplot(2,2,1)
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides making stacking interaction');
subplot(2,2,2)
boxplot(Response(i));
subplot(2,2,3)
i = find(Data(:,c) == 0);
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides not making a stacking interaction');
subplot(2,2,4)
boxplot(Response(i));
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf



c = 6;
i = find(Data(:,c) > 0);
subplot(2,2,1)
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides making base-phosphate interaction');
subplot(2,2,2)
boxplot(Response(i));
subplot(2,2,3)
i = find(Data(:,c) == 0);
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides not making a base-phosphate interaction');
subplot(2,2,4)
boxplot(Response(i));
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf



c = 7;
i = find(Data(:,c) > 0);
subplot(2,2,1)
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides near proteins');
subplot(2,2,2)
boxplot(Response(i));
subplot(2,2,3)
i = find(Data(:,c) == 0);
hist(Response(i),Range)
title('Nucleotide conservation - nucleotides not near proteins');
subplot(2,2,4)
boxplot(Response(i));
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,Data(:,4));
ylabel(Resp);
title('Number of cWW basepairs made')
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,Data(:,5));
ylabel(Resp);
title('Number of non-cWW basepairs made')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,Data(:,6));
ylabel(Resp);
title('Number of base-phosphate interactions made')
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,Data(:,8));
ylabel(Resp);
title('Number of stacking interactions made')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4) Data(:,5)});
ylabel(Resp);
title('Number of cWW basepairs above non-cWW basepairs')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4) Data(:,6)});
ylabel(Resp);
title('Number of cWW basepairs above base-phosphates')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,5) Data(:,6)});
ylabel(Resp);
title('Number of non-cWW basepairs above base-phosphates')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4) Data(:,8)});
ylabel(Resp);
title('Number of cWW basepairs above number of stacks')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,8) Data(:,6)});
ylabel(Resp);
title('Number of stacks above number of base-phosphates')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4)>0 Data(:,5)>0 Data(:,6)>0});
ylabel(Resp);
title('Watson-Crick versus non-Watson-Crick versus BPh')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,8)>0 Data(:,6)>0});
ylabel(Resp);
title('stack versus base-phosphates')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,8)>0 Data(:,6)});
ylabel(Resp);
title('stack versus number of base-phosphates')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4)>0 Data(:,5)>0 Data(:,7)>0});
ylabel(Resp);
title('cWW versus non-cWW versus near protein')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
boxplot(Response,{Data(:,4)>0 Data(:,5)>0 Data(:,8)>0});
ylabel(Resp);
title('cWW versus non-cWW versus stacking')
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

figure(fn)
clf
scatter(Data(:,3),Data(:,11),6,'filled')
title('Graph of entropy versus conservation percentage');
xlabel('Conservation percentage');
ylabel('Entropy');
Num = num2str(fn);
fn = fn + 1;
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num '.pdf'],'pdf');
saveas(gcf,['Conservation' filesep 'ConservationPrediction-' Num 'LowRes.png'],'png');

% break



% Standardize the predictors if you want to see how strong each type of
% interaction contributes (regardless of the number of such interactions

[n,v] = size(Data(:,Variables));
m = mean(Data(:,Variables));                     % mean
s = std(Data(:,Variables));                      % standard deviation
m = 0*m;


StData = (Data(:,Variables) - ones(n,1)*m) ./ (ones(n,1)*s); % standardized
StData(:,v) = ones(n,1);                         % constant column


fprintf('Linear regression model for all nucleotides:\n');
fprintf('%s = \n', Resp);
fprintf('%7.4f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.4f * %7.4f * (%15s - %7.4f)\n', abs(B(v)), s(v), VarName{Variables(v)},m(v));
end

Variables = Variables(2:end);
StData = StData(:,2:end);
s = s(2:end);
m = m(2:end);

i = find(Data(:,4) > 0);             % cWW pairs
[B,BINT,R,RINT,STATS] = regress(Response(i),StData(i,:));
fprintf('Linear regression model for nucleotides making a cWW pair:\n');
fprintf('%s = \n', Resp);
fprintf('%7.4f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.4f * %7.4f * (%15s - %7.4f)\n', abs(B(v)), s(v), VarName{Variables(v)},m(v));
end

i = find(Data(:,4) == 0);             % not making a cWW pair
[B,BINT,R,RINT,STATS] = regress(Response(i),StData(i,:));
fprintf('Linear regression model for nucleotides not making a cWW pair:\n');
fprintf('%s = \n', Resp);
fprintf('%7.4f (constant term)\n', B(end));
for v = 1:(length(Variables)-1),
  if B(v) > 0,
    fprintf('      + ');
  else
    fprintf('      - ');
  end
  fprintf('%7.4f * %7.4f * (%15s - %7.4f)\n', abs(B(v)), s(v), VarName{Variables(v)},m(v));
end

break

% ---------------------------------------- study conservation of phosphate acc

i = find(Data(:,12) > 0);                % phosphate acceptors
D = Data(i,:);

D = sortrows(D,[4 3]);

for i = 1:length(D(:,1)),
  f = D(i,1);
  fprintf('%8.4f ', D(i,3));
  zListNucleotideInteractions(File(f),D(i,2));

  if D(i,4) > 0,                         % makes a cWW pair
    j = find(abs(fix(File(f).Edge(D(i,2),:))) == 1);
    if ~isempty(j),
      fprintf('         ');
      zListNucleotideInteractions(File(f),j);
    end
  end

  e = abs(fix(File(f).Edge(D(i,2),:)));
  j = find( (e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    fprintf('         Conservation of stacking partners %8.4f\n', mean(c));
  end
end

% ---------------------------------------- study conservation of cWW pairs

figure(1)
clf
Range    = 1+[0:2:98];

subplot(2,2,1)
i = find(Data(:,4) > 0);                 % makes a cWW pair
hist(Data(i,3),Range);
xlabel('Conservation percentage');
title('All cWW pairs');

subplot(2,2,2)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0));  % makes a cWW pair
hist(Data(i,3),Range);
xlabel('Conservation percentage');
title('cWW pair but no non-cWW');

subplot(2,2,3)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,12) > 0));  % makes a cWW pair
hist(Data(i,3),Range);
xlabel('Conservation percentage');
title('cWW, no non-cWW, acceptor');

subplot(2,2,4)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,10) > 0));  % makes a cWW pair
hist(Data(i,3),Range);
xlabel('Conservation percentage');
title('cWW, no non-cWW, donor');

% ---------------------------------------- study conservation of cWW pairs

figure(2)
clf
Range    = 1:4;

subplot(2,2,1)
i = find(Data(:,4) > 0);                 % makes a cWW pair
hist(Data(i,17),Range);
xlabel('Base A C G U');
title('All cWW pairs');

subplot(2,2,2)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0));  % makes a cWW pair
hist(Data(i,17),Range);
xlabel('Base A C G U');
title('cWW pair but no non-cWW');

subplot(2,2,3)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,12) > 0));  % makes a cWW pair
hist(Data(i,17),Range);
xlabel('Base A C G U');
title('cWW, no non-cWW, acceptor');

subplot(2,2,4)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,10) > 0));  % makes a cWW pair
hist(Data(i,17),Range);
xlabel('Base A C G U');
title('cWW, no non-cWW, donor');

% ---------------------------------------- study conservation of cWW pairs

figure(3)
clf
Range    = 1+[0:2:98];

subplot(2,2,1)
i = find(Data(:,4) > 0);                 % makes a cWW pair
for a = 1:length(i),
  f = Data(i(a),1);
  e = abs(fix(File(f).Edge(Data(i(a),2),:)));
  j = find( (e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    m(a) = mean(c);
  else
    m(a) = 0;
  end
end  

subplot(2,2,2)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0));  % makes a cWW pair
clear m
for a = 1:length(i),
  f = Data(i(a),1);
  e = abs(fix(File(f).Edge(Data(i(a),2),:)));
  j = find( (e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    m(a) = mean(c);
  else
    m(a) = 0;
  end
end  
plot(Data(i,3),m,'.')
xlabel('Conservation of NT');
ylabel('Conservation of stacking partners');
title('cWW pair but no non-cWW');

subplot(2,2,3)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,12) > 0));  % makes a cWW pair
clear m
for a = 1:length(i),
  f = Data(i(a),1);
  e = abs(fix(File(f).Edge(Data(i(a),2),:)));
  j = find( (e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    m(a) = mean(c);
  else
    m(a) = 0;
  end
end  
plot(Data(i,3),m,'.')
xlabel('Conservation of NT');
ylabel('Conservation of stacking partners');
title('cWW, no non-cWW, acceptor');

subplot(2,2,4)
i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,10) > 0));  % makes a cWW pair
clear m
for a = 1:length(i),
  f = Data(i(a),1);
  e = abs(fix(File(f).Edge(Data(i(a),2),:)));
  j = find( (e > 19) .* (e < 25));                % stacking partners
  c = File(f).Conservation(j);
  c = c(find(c > 0));
  if ~isempty(c),
    m(a) = mean(c);
  else
    m(a) = 0;
  end
end  
plot(Data(i,3),m,'.')
xlabel('Conservation of NT');
ylabel('Conservation of stacking partners');
title('cWW, no non-cWW, donor');


% ---------------------------------------- study conservation of cWW pairs

figure(3)
clf
Range    = 1+[0:2:98];

for v = 1:4,
  subplot(2,2,v);
  switch v,
  case 1, 
    i = find(Data(:,4) > 0);                 % makes a cWW pair
  case 2,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0));  % makes a cWW pair
  case 3,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,12) > 0));  % makes a cWW pair
  case 4,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,10) > 0));  % makes a cWW pair
  end

  clear m
  for a = 1:length(i),
    f = Data(i(a),1);
    e = abs(fix(File(f).Edge(Data(i(a),2),:)));
    j = find( (e > 19) .* (e < 25));                % stacking partners
    c = File(f).Conservation(j);
    c = c(find(c > 0));
    if ~isempty(c),
      m(a) = mean(c);
    else
      m(a) = 0;
    end
  end  

  plot(Data(i,3),m,'.')
  xlabel('Conservation of NT');
  ylabel('Conservation of stacking partners');

  switch v,
  case 1, title('All cWW pairs');
  case 2, title('cWW pair but no non-cWW');
  case 3, title('cWW, no non-cWW, acceptor');
  case 4, title('cWW, no non-cWW, donor');
  end
end

% ---------------------------------------- study conservation of cWW pairs

figure(4)
clf
Range    = 1+[0:2:98];

for v = 1:4,
  subplot(2,2,v);
  switch v,
  case 1, 
    i = find(Data(:,4) > 0);                 % makes a cWW pair
  case 2,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0));  % makes a cWW pair
  case 3,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,12) > 0));  % makes a cWW pair
  case 4,
    i = find((Data(:,4) > 0) .* (Data(:,5) == 0) .* (Data(:,10) > 0));  % makes a cWW pair
  end

  clear m
  for a = 1:length(i),
    f = Data(i(a),1);
    e = abs(fix(File(f).Edge(Data(i(a),2),:)));
    j = find( (e > 19) .* (e < 25));                % stacking partners
    c = File(f).Conservation(j);
    c = c(find(c > 0));
    if ~isempty(c),
      m(a) = mean(c);
    else
      m(a) = 0;
    end
  end  
  hist(m,Range);
  xlabel('Conservation of stacking partners');

  mean(m)

  switch v,
  case 1, title('All cWW pairs');
  case 2, title('cWW pair but no non-cWW');
  case 3, title('cWW, no non-cWW, acceptor');
  case 4, title('cWW, no non-cWW, donor');
  end
end

figure(5)
clf
plot(Data(:,3))
ylabel('Conservation percentage')
xlabel('Position along the chain')

figure(6)
clf
plot(Data(:,3),Data(:,18),'.')
ylabel('Conservation percentage')
xlabel('Position along the chain')

