% This program needs revisions!  This version has garbage from a specific problem!

% pGroupSequences clusters sequences and displays groups together

function [Group,Group2] = pGroupSequences(Node,Sequence,NumGroups,Group,Group2)

if nargin < 3,
  NumGroups = min(ceil(length(Sequence)/2),8);
end

if nargin < 5,

  Dist = zeros(length(Sequence));

  for k=1:length(Sequence),
    a = diff([cat(2,Sequence(k).TraceInfo(:).mp) 0]);
    for j=1:k-1,
      b = diff([cat(2,Sequence(j).TraceInfo(:).mp) 0]);
      Dist(j,k) = sqrt(sum((a-b).^2));
      Dist(k,j) = Dist(j,k);
    end
  end

  [Group,Group2] = pCluster(Dist);
else

[y,i] = sort(Group2(NumGroups,:));

diary(['zParB_' num2str(NumGroups) '_group.txt']);

fprintf('5S RNA nucleotides 71 to 110 with %4d groups\n\n', NumGroups);
fprintf('Aligned according to bacterial model\n');

zParseAlignment(Node,Sequence(i),Group(NumGroups,i));

fprintf('\n\n5S RNA nucleotides 71 to 110 with %4d groups\n', NumGroups);
fprintf('----------- Each group is aligned separately ------------\n');

for g = 1:NumGroups,
  fprintf('\n');
  j = find(Group(NumGroups,i) == g);
  zParseAlignment(Node,Sequence(i(j)),Group(NumGroups,i(j)));
end

diary off

end