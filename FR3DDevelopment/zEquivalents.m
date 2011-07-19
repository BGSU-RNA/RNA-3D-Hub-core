% zEquivalents identifies equivalent 3D structures and sorts them by some criterion

function [t,n] = zEquivalents(t,n,Criterion,Verbose)

% load PDBInfo

if nargin < 3,
  Criterion = 5;
end

if nargin < 4,
  Verbose = 1;
end

[a,b] = size(t);                          % a is the number of files

[y,i] = sort(-n(:,4));                    % decreasing order of size
t = t(i,:);                               % re-order files
n = n(i,:);

for i = 1:a,
  if isempty(t{i,10}),
    t{i,10} = t{i,1};                     % these must be in a group alone
  end
end

U = unique(t(:,10));                      % names of groups

Order = [(1:a)' zeros(a,1)];              % nominal order of files

for u = 1:length(U),                      % loop through the groups
  j = find(ismember(t(:,10),U{u}));       % find the others in the same group
  i = min(j);                             % use first member of group as ID
  for jj = 1:length(j),                   % work through members of the group
    Order(j(jj),1) = i;                   % rank group all the same at first
    Order(j(jj),2) = datenum(t{j(jj),4}, 'yyyy-mm-dd');
    Order(j(jj),3) = n(j(jj),1);          % resolution
    Order(j(jj),4) = -n(j(jj),2);         % number of nucleotides
    Order(j(jj),5) = -n(j(jj),3);         % number of pairs
    if n(j(jj),2) > 0,
      Order(j(jj),6) = n(j(jj),3)/n(j(jj),2);    
                                          % pairs per nucleotide
    else
      Order(j(jj),6) = 0.1;    
    end
    if n(j(jj),2) > 0,
      Order(j(jj),7) = n(j(jj),8)/n(j(jj),7);    
                                          % pairs per nucleotide
    else
      Order(j(jj),7) = 0.1;    
    end
  end
end

switch mod(Criterion,10)
  case 1, [y,p] = sortrows(Order,[1 -6 3 -4 -5]);
  case 2, [y,p] = sortrows(Order,[1 3 -6]);
  case 3, [y,p] = sortrows(Order,[1 -4 3]);
  case 4, [y,p] = sortrows(Order,[1 -5 -6]);
  case 5, [y,p] = sortrows(Order,[1 -6 3 -2]);
  case 5, [y,p] = sortrows(Order,[1 -6 3 -2]);
  case 6, [y,p] = sortrows(Order,[1 -7 3 -2]);
end

t = t(p,:);
n = n(p,:);
Order = Order(p,:);

if Verbose > 0,
  for i = 1:a,
    fprintf('%4s %4s %4d %4d %4d %10.6f %10.6f\n', t{i,1}, t{i,10}, n(i,2), n(i,4), Order(i,1), Order(i,6), Order(i,3));
  end
end
