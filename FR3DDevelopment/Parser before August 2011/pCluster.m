
function [Group,Group2] = zCluster(Dist)

L = length(Dist);

Max = max(max(Dist));

ODist = Dist;                             % original distances

NG = L;                                   % number of groups right now

Group(NG,:) = 1:L;                        % group numbers with N groups

while NG > 1,
  [N,a] = min(Dist+eye(length(Dist))*Max);% find smallest non-zero distance
  [M,b] = min(N);

  c = a(b);                               % (c,d) is the minimal pair
  d = b;

  if Group(NG,c) ~= Group(NG,d),          % items come from different groups
    g = Group(NG,c);
    Group(NG-1,:) = Group(NG,:);          % keep other group numbers
    Group(NG-1,find(Group(NG,:) == Group(NG,d))) = g;
    NG = NG - 1;                          % number of groups is smaller
fprintf('Number of groups %4d\n', NG);
drawnow
    A = zeros(1,L);                        % now renumber groups
    for j=1:L,
      A(Group(NG,j)) = A(Group(NG,j)) + 1;  % count number in each group
    end
    [y,i] = sort(-A);                       % sort in decreasing order
    k(i) = 1:L;                            % k is inverse of i
    for j=1:L,
      Group(NG,j) = k(Group(NG,j));      
    end
  end

  Dist(c,d) = Max;                        % do not consider this pair again
  Dist(d,c) = Max;
end

% add distances to other members of the group, for subclassification

for NG = 1:L,
  for j=1:L,
    p = find(Group(NG,:) == Group(NG,j));
    Group2(NG,j) = Group(NG,j) + sum(ODist(j,p))/(100*Max*L) + Group(min(NG+10,L),j)/(2*Max*L);
  end
end