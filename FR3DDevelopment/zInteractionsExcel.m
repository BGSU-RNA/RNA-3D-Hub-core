
for f = 1:length(File),

clear T
E = fix(File(f).Edge);
BPh = File(f).BasePhosphate;

W = [1:6 -1 -2];
H = [7 8 9 10 -3 -4 -7 -8];
S = [11 12 -5 -6 -9 -10 -11 -12];

WBPh = [1 2 5 11 12 13 19 15];
HBPh = [3 4 6 7 8 9 18 14 16 17];
SBPh = [1 10];

for i = 1:length(File(f).NT),
  N = File(f).NT(i);
  T{i,1} = N.Base;
  T{i,2} = N.Number;
  for k = 3:14,
    T{i,k} = '';
  end
  
  j = find(E(i,:) );
  for k = 1:length(j),
    if any(E(i,j(k)) == W),
      T{i,3} = [T{i,3} zEdgeText(E(i,j(k)),0,N.Code,File(f).NT(j(k)).Code)];
      T{i,4} = [T{i,4} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
    if any(E(i,j(k)) == H),
      T{i,5} = [T{i,5} zEdgeText(E(i,j(k)),0,N.Code,File(f).NT(j(k)).Code)];
      T{i,6} = [T{i,6} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
    if any(E(i,j(k)) == S),
      T{i,7} = [T{i,7} zEdgeText(E(i,j(k)),0,N.Code,File(f).NT(j(k)).Code)];
      T{i,8} = [T{i,8} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
  end

  BPh(i,i) = 0;                           % eliminate self interactions

  j = find(BPh(i,:) );
  for k = 1:length(j),
    if any(BPh(i,j(k)) == WBPh),
      T{i,9} = [T{i,9} zBasePhosphateText(BPh(i,j(k)))];
      T{i,10} = [T{i,10} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
    if any(BPh(i,j(k)) == HBPh),
      T{i,11} = [T{i,11} zBasePhosphateText(BPh(i,j(k)))];
      T{i,12} = [T{i,12} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
    if any(BPh(i,j(k)) == SBPh),
      T{i,13} = [T{i,13} zBasePhosphateText(BPh(i,j(k)))];
      T{i,14} = [T{i,14} File(f).NT(j(k)).Base File(f).NT(j(k)).Number ' '];
    end
  end

end

clear H

H{1} = 'Base';
H{2} = 'Number';
H{3} = 'W edge with';
H{4} = '';
H{5} = 'H edge with';
H{6} = '';
H{7} = 'S edge with';
H{8} = '';
H{9} = 'W edge with';
H{10} = '';
H{11} = 'H edge with';
H{12} = '';
H{13} = 'S edge with';
H{14} = '';

T = [H; T];
T = T(:,[1 2 3 4 9 10 5 6 11 12 7 8 13 14]);

xlswrite([File(f).Filename '_InteractionList.xls'],T);

end
