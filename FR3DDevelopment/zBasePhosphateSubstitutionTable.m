
BPhList = [1 2 3 4 5 6 7 18 8 9 10 11 12 19 13 14 15 16 17];

for b = 1:17,
%  BPh = BPhList(b);
  BPh = b;
  [Distance, Score, MinBPh] = zBasePhosphateGeometry(BPh);
  [T,B] = zBasePhosphateText(BPh);
  fprintf('%s%s nearest distances: ', B, T);
  for c = 1:4,
    [T,B] = zBasePhosphateText(MinBPh(c));
    fprintf('%s%s %8.4f    ', B, T, Distance(c));
  end
  fprintf('\n');
end

for b = 1:17,
%  BPh = BPhList(b);
  BPh = b;
  [Distance, Score, MinBPh] = zBasePhosphateGeometry(BPh);
  [T,B] = zBasePhosphateText(BPh);
  fprintf('%s%s substitution scores: ', B, T);
  for c = 1:4,
    [T,B] = zBasePhosphateText(MinBPh(c));
    fprintf('%s %8.4f    ', B, Score(c));
  end
  fprintf('\n');
end



