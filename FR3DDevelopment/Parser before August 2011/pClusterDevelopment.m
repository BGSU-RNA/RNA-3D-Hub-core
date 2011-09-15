
L = rand(10,2);

for i=1:(length(L)-1)
  for j=(i+1):length(L) 
    Dist(i,j)=sum((L(i,:)-L(j,:)).^2);       % L^2 distance between them
    Dist(j,i)=Dist(i,j);
  end
end

n1 = length(Dist);

Max = max(max(Dist));
Sum=sum(Dist')';

NG = n1;                                   % number of groups right now

Group(NG,:) = 1:n1;                        % group numbers with N groups

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

    A = zeros(1,n1);                        % now renumber groups
    for j=1:n1,
      A(Group(NG,j)) = A(Group(NG,j)) + 1;  % count number in each group
    end
    [y,i] = sort(-A);                       % sort in decreasing order
    k(i) = 1:n1;                            % k is inverse of i
    for j=1:n1,
      Group(NG,j) = k(Group(NG,j));      
    end
  end

  Dist(c,d) = Max;                        % do not consider this pair again
  Dist(d,c) = Max;
end

figure(1)

for i=n1:-1:1,
  fprintf('%2d groups\n',i);
  fprintf('%2d ', Group(i,:));
  fprintf('\n');
  scatter(L(:,1),L(:,2),5*ones(n1,1),Group(i,:),'filled');
  pause
end

