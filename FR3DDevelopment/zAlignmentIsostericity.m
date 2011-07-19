% zAlignmentIsostericity(File,i1,i2) looks at each basepair that can be inferred in File(2) from the alignment i1, i2 and the basepairs in File(1) and calculates the IDI values between the exemplar of the basepair found in File(1) with the exemplar of the base combination from File(2).  It does not use the actual basepair found in File(2), only the base combination that is aligned to the basepair in File(1).

function [IDI] = zAlignmentIsostericity(File,i1,i2)

load PairExemplars

onetotwo = zeros(1,length(File(1).NT));

onetotwo(i1) = i2;

E = File(1).Edge;

[a,b,c] = find(triu(E) .* (abs(E) > 0) .* (abs(E) < 13));  % pairs of basepairs

IDI = [];

for k = 1:length(a),
  x = onetotwo(a(k));
  y = onetotwo(b(k));
  if x > 0 && y > 0,   % both indices are aligned
% [E(a(k),b(k)) File(1).NT(a(k)).Code File(1).NT(b(k)).Code]
    idimatrix = zExemplarIDI(E(a(k),b(k)),File(1).NT(a(k)).Code,File(1).NT(b(k)).Code,ExemplarIDI);
    h = idimatrix(File(2).NT(x).Code,File(2).NT(y).Code);
    if isnan(h),
      h = 20;
    end

    if h > 5,
      fprintf('Large IDI (%7.4f) with %s%4s - %s%4s %s aligned to %s%4s - %s%4s %s.\n', h, ...
              File(2).NT(x).Base,File(2).NT(x).Number,File(2).NT(y).Base,File(2).NT(y).Number, zEdgeText(File(2).Edge(x,y)), ...
              File(1).NT(a(k)).Base,File(1).NT(a(k)).Number,File(1).NT(b(k)).Base,File(1).NT(b(k)).Number, zEdgeText(File(1).Edge(a(k),b(k))));
    end

    IDI = [IDI; h];
  end
end

hist(IDI,50)
