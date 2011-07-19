
% loading code at the end

clear d

for cc = 7:7,

fprintf('Analyzing the %s family\n', zEdgeText(cc));

d{1} = Seq_Alignment_Extraction_23S{cc};
d{2} = Seq_Alignment_Extraction_16S{cc};
d{3} = Seq_Alignment_Extraction_5S{cc};

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

Pairs = {'AA','CA','GA','UA','AC','CC','GC','UC','AG','CG','GG','UG','AU','CU','GU','UU'};



S = 0;
T = 0;

for i = 1:length(d),
  [s,t] = size(d{i});
  S = max(s,S);
  T = max(t,T);
end

Counts = zeros(16,3*T);      % a 16 x 3T matrix of counts down columns

Structure = {'23S','16S','5S'};

for i = 1:length(d),
  fprintf('Accumulating counts from %s\n', Structure{i});
  [s,t] = size(d{i});
  for k = 1:t,                  % go through the columns
    for j = 1:s,                  % go down each row
      m = find(ismember(Pairs,d{i}{j,k}));
      if ~isempty(m),
        Counts(m,(i-1)*T+k) = Counts(m,(i-1)*T+k) + 1;
      end
    end
  end
end

% ------------------------------- Symmetrize for cWW, tWW, ... families

NumPairs = 16;
pc = 1:16;

if any(cc == [1 2 7 8]),

Counts(5,:)  = Counts(5,:)  + Counts(2,:);  % AC
Counts(9,:)  = Counts(9,:)  + Counts(3,:);  % AG
Counts(10,:) = Counts(10,:) + Counts(7,:);  % GC
Counts(13,:) = Counts(13,:) + Counts(4,:);  % AU
Counts(14,:) = Counts(14,:) + Counts(8,:);  % CU
Counts(15,:) = Counts(15,:) + Counts(12,:); % GU
Counts(2,:)  = 0 * Counts(2,:);
Counts(3,:)  = 0 * Counts(3,:);
Counts(7,:)  = 0 * Counts(7,:);
Counts(4,:)  = 0 * Counts(4,:);
Counts(8,:)  = 0 * Counts(8,:);
Counts(12,:) = 0 * Counts(12,:);

NumPairs = 10;
pc = [1 5 6 9 10 11 13 14 15 16];

end

% ------------------------------- Remove empty columns (due to sizes of 5S,16S)

Total = sum(Counts,1);

i = find(Total);
T = length(i);

fprintf('Found %d %s basepairs in the conserved core\n', T, zEdgeText(cc));

Counts = Counts(:,i);
Total  = Total(i);

% ------------------------------- Compute estimated percentages from each col

Perc = Counts ./ (ones(16,1) * Total);
Numb = sum(Counts,2);             % number of observations of each type

% ------------------------------- Save for later

BPCountsFromSequences{cc} = Counts;

% ------------------------------- Estimates using all columns equally

pp = sum(Perc,2) / T;           % weigh each column equally to est percentages
for m = 1:16,
  if sum(Counts(m,:)) > 0,
    fprintf('%2s %6d %7.4f\n', Pairs{m}, Numb(m), 100*pp(m));
  end
end


if T > 10,

% ------------------------------- Bootstrap estimates of probabilities

B = 1000000;                     % number of bootstrap estimates

fprintf('Conducting %d bootstrap calculations of the proportions,\n', B);

p = zeros(16,B);                % to store each estimated proportion

for b = 1:B,
  c = ceil(T * rand(1,T));      % choose T columns randomly, iid
  p(:,b) = sum(Perc(:,c),2) / T;
end

% ------------------------------ Choose intervals to get good family rate

FamSuccess = 0;
a = 0.000000001;               % lower limit on q
b = 0.05;                      % upper limit on q
q = (a+b)/2;
n = 1;

while ((FamSuccess < 0.95)||(FamSuccess > 0.9505)) && (n < 50),
  l = quantile(p,  q,2);       % lower limit of 95% confidence interval
  u = quantile(p,1-q,2);       % upper limit of 95% confidence interval

  % ------------------------------ Calculate the family-wise error rate

  OKInt = (p >= l*ones(1,B)) .* (p <= u*ones(1,B));

  AllOK = (sum(OKInt,1) == 16);
  FamSuccess = sum(AllOK)/B;

  fprintf('Level %16.14f AA left %8.6f AA right %8.6f Family success %8.6f\n', q, l(1), u(1), FamSuccess);

  if FamSuccess < 0.95, 
    b = q;
  elseif FamSuccess > 0.95,
    a = q;
  end

  q = (a+b)/2;
  n = n + 1;
end

fprintf('After %d bootstrap calculations of the proportions,\n', B);

for m = 1:16,
  if sum(Counts(m,:)) > 0,
    fprintf('%2s estimate %7.4f. Confidence interval (%7.4f, %7.4f)\n', Pairs{m}, 100*pp(m), 100*l(m), 100*u(m));
  end
end

CI{cc} = [Numb 100*pp 100*l 100*u];

else

% Use the Quesenberry formula to calculate

L = gaminv(1-0.05/NumPairs,1/2,2);    % chi-squared cutoff from Genz
N = T;                                % number of observations
n = pp*N;                             % equivalent number of instances

l = (L + 2*n - sqrt(L*(L + 4*n.*(N-n)/N)))/(2*(L+N));
u = (L + 2*n + sqrt(L*(L + 4*n.*(N-n)/N)))/(2*(L+N));

for m = 1:16,
  if sum(Counts(m,:)) > 0,
    fprintf('%2s estimate %7.4f. Confidence interval (%7.4f, %7.4f)\n', Pairs{m}, 100*pp(m), 100*l(m), 100*u(m));
  end
end

CI{cc} = [Numb 100*pp 100*l 100*u];

end

end

break

load 23S_basepairs_from_Seqs.mat
load 16S_basepairs_from_Seqs.mat
load 5S_basepairs_from_Seqs.mat

break

% Code to overwrite cells in an Excel table

Lett = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% ---------------------------------------------- Write 16x16 tables

ExcelTemp

