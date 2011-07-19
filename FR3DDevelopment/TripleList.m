
c = [1 2 3 -3 4 -4 5 -5 6 -6 7 8 9 -9 10 -10 11 12];
c = [1:12 -3 -4 -5 -6 -9 -10];
h = [];

for i = c,
  for j = c,
    e = zEdgeText(i);
    f = zEdgeText(j);

    if lower(e(3)) ~= lower(f(2)),
      h = [h; [e(1:3) f(1:3)]];
%      fprintf('%s%s\n',e(1:3), f(1:3));
    end
  end
end

R = h(:,[4 6 5 1 3 2]);

OK = ones(length(h(:,1)),1);

for i = length(h(:,1)):-1:1,
  for j = 1:i,
    if all(h(i,:) == R(j,:)),
      OK(i) = 0;
    end
  end
end

h = h(find(OK),:)

fprintf('Found %d different triples\n', length(find(OK)));  