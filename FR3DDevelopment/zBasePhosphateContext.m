
File = zAddNTData('2aw4');

if 1 < 0,
  File = zAddNTData('1s72');

  File = xAnnotateWithKnownMotifs(File,1,0);
  File = xAnnotateWithKnownMotifs(File,1,0,{'BPh_pair_triple.mat'});
end

S.File = File;
S.Query.Geometric = 0;

fprintf('Base-phosphate report for %s.\n', File.Filename);

BPh = File.BasePhosphate;

for m = 1:length(BPh(1,:)),
  BPh(m,m) = 0;                             % omit self interactions
end

[i,j,k] = find(BPh .* (BPh > 0) .* (BPh < 100));

fprintf('Found %d distinct base-phosphate interactions.\n', length(i));
fprintf('Found %d distinct bases and %d distinct phophates.\n', length(unique(i)), length(unique(j)));

figure(1)
clf
d = diff(sort(unique(i)));
hist(d,1:max(d));
hold on
p = length(unique(i)) / File.NumNT;            % percentage making BPh
plot(1:max(d),length(unique(i))*p*(1-p).^(1:max(d)));

figure(2)
clf
d = diff(sort(unique(j)));
hist(d,1:max(d));
hold on
p = length(unique(j)) / File.NumNT;            % percentage making BPh
plot(1:max(d),length(unique(j))*p*(1-p).^(1:max(d)));

figure(3)
clf
d = diff(sort(unique([i j])));
hist(d,1:max(d));
hold on
p = length(unique([i j])) / File.NumNT;            % percentage making BPh
plot(1:max(d),length(unique([i j]))*p*(1-p).^(1:max(d)));

IR = [];
Edge = [];
for m = 1:length(i),
  IR(m,1) = File.Range(i(m),j(m));
  Edge(m,1) = File.Edge(i(m),j(m));
end

BPCat = [2 6 7 0 6 7 8 9 0 1 3 4 5 0 5 9 0 8 4]';  % updated 8-19-2008

Data = [i j k IR i-j Edge min(5,IR) min(abs(i-j),3) fix(mod(Edge,100)) BPCat(k)];

% columns:
%  1 index of the base
%  2 index of the phosphate
%  3 code of the base-phosphate interaction
%  4 interaction range
%  5 signed distance in sequence
%  6 pairwise interaction between the two
%  7 IR limited to size 5
%  8 distance in sequence limited to 3
%  9 edge, not distinguishing near and true
% 10 BPh category 0 to 9

% ---------------- adjacent and stacked

r = find((Data(:,8) == 1).*(abs(Data(:,6)) > 19).*(abs(Data(:,6))<24));
%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

fprintf('Of the %d interactions, %3d (%6.2f%%) are between nucleotides that are adjacent in the nucleotide sequence and are stacked on one another.\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% ---------------- adjacent and not stacked

r = find(Data(:,8) == 1);
%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

fprintf('Of the %d interactions, %3d (%6.2f%%) are between nucleotides that are adjacent in the nucleotide sequence and are not stacked on one another.\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% ---------------- both are part of a GNRA hairpin
Te = 'HL=';
r = [];
for a = 1:length(Data(:,1)),
  if ~isempty(File.Nucl(Data(a,1)).Motif) && ~isempty(File.Nucl(Data(a,2)).Motif),
    if ~isempty(strfind(File.Nucl(Data(a,1)).Motif(1).Name,Te)) && ~isempty(strfind(File.Nucl(Data(a,2)).Motif(1).Name,Te)),
      r = [r; a];
    end
  end
end

fprintf('Of the %d interactions, %3d (%6.2f%%) are between nucleotides in a GNRA, T-loop, or UNCG hairpin.\n', length(i), length(r), 100*length(r)/length(i));

%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% ---------------- BPh_pair_triple

Te = 'BPh_pair_triple';

r = [];
for a = 1:length(Data(:,1)),
  if ~isempty(File.Nucl(Data(a,1)).Motif) && ~isempty(File.Nucl(Data(a,2)).Motif),
    if ~isempty(strfind(File.Nucl(Data(a,1)).Motif(1).Name,Te)) && ~isempty(strfind(File.Nucl(Data(a,2)).Motif(1).Name,Te)),
      r = [r; a];
    end
  end
end

%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

fprintf('Of the %d interactions, %3d (%6.2f%%) are part of a nested triple in which X makes a BPh with Y, X makes a nested basepair with Z, and Y and Z are sequential.\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% ---------------- Internal loop

Te = 'IL_';

r = [];
for a = 1:length(Data(:,1)),
  if ~isempty(File.Nucl(Data(a,1)).Motif) && ~isempty(File.Nucl(Data(a,2)).Motif),
    if ~isempty(strfind(File.Nucl(Data(a,1)).Motif(1).Name,Te)) && ~isempty(strfind(File.Nucl(Data(a,2)).Motif(1).Name,Te)),
      r = [r; a];
    end
  end
end

%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

fprintf('Of the %d interactions, %3d (%6.2f%%) are part of an internal loop such as sarcin or a kink turn.\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% --------------------

r = find((Data(:,4) == 0) .* (abs(Data(:,6) == 8)));
%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

fprintf('Of the %d interactions, %3d (%6.2f%%) are between nucleotides that are also making a tHH basepair within a single stem loop (nested interaction).\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

% -------------------- Other nested interactions

r = find(Data(:,4) == 0);
%xDisplayCandidates(File,[Data(r,[1 2]) ones(length(r),1)])

for a = 1:length(r),
  fprintf('%4d Base %4s makes interactions ', a, File.NT(Data(a,1)).Number);
  fprintf('%6s with %4s, ', zBasePhosphateText(Data(a,3)), File.NT(Data(a,2)).Number);
  E = abs(File.Edge(Data(a,1),:));
  q = find((E > 0) .* (E < 30));
  [y,i] = sort(abs(E(q)));
  q = q(i);
  for b = 1:length(q),
    fprintf('%5s with %4s, ', zEdgeText(File.Edge(Data(a,1),q(b))), File.NT(q(b)).Number);
  end
  fprintf('\n');

  fprintf('%4d Phos %4s makes interactions ', a, File.NT(Data(a,2)).Number);
  fprintf('%6s with %4s, ', zBasePhosphateText(-Data(a,3)), File.NT(Data(a,1)).Number);
  E = abs(File.Edge(Data(a,2),:));
  q = find((E > 0) .* (E < 30));
  [y,i] = sort(abs(E(q)));
  q = q(i);
  for b = 1:length(q),
    fprintf('%5s with %4s, ', zEdgeText(File.Edge(Data(a,2),q(b))), File.NT(q(b)).Number);
  end
  fprintf('\n');

end




fprintf('Of the %d interactions, %3d (%6.2f%%) are between nucleotides that are interacting within a single stem loop (nested interaction).\n', length(i), length(r), 100*length(r)/length(i));

g = setdiff(1:length(Data(:,1)),r);
Data = Data(g,:);                        % keep only the others

r = 1:length(Data(:,1));

fprintf('Of the %d interactions, %3d (%6.2f%%) remain to be classified.\n', length(i), length(r), 100*length(r)/length(i));


for a = 1:length(r),
  fprintf('%4d Base %4s makes interactions ', a, File.NT(Data(a,1)).Number);
  fprintf('%6s with %4s, ', zBasePhosphateText(Data(a,3)), File.NT(Data(a,2)).Number);
  E = abs(File.Edge(Data(a,1),:));
  q = find((E > 0) .* (E < 30));
  [y,i] = sort(abs(E(q)));
  q = q(i);
  for b = 1:length(q),
    fprintf('%5s with %4s, ', zEdgeText(File.Edge(Data(a,1),q(b))), File.NT(q(b)).Number);
  end
  fprintf('\n');

  fprintf('%4d Phos %4s makes interactions ', a, File.NT(Data(a,2)).Number);
  fprintf('%6s with %4s, ', zBasePhosphateText(-Data(a,3)), File.NT(Data(a,1)).Number);
  E = abs(File.Edge(Data(a,2),:));
  q = find((E > 0) .* (E < 30));
  [y,i] = sort(abs(E(q)));
  q = q(i);
  for b = 1:length(q),
    fprintf('%5s with %4s, ', zEdgeText(File.Edge(Data(a,2),q(b))), File.NT(q(b)).Number);
  end
  fprintf('\n');

end







break

fprintf('Of these, %d have interaction range 0, %d have interaction range 1 to 5, and %d have interaction range above 5.\n', length(find(IR == 0)), length(find((IR > 0) .* (IR < 6))), length(find(IR > 5)));


Data = sortrows(Data,[7 9 8 10]);

S.Candidates = [Data(:,[1 2]) ones(length(i),1)];
S.Discrepancy = ones(length(i),1);
xListCandidates(S);

xDisplayCandidates(File,S.Candidates);

break






fprintf('Pairs having interaction range 0:\n');

r = find(Data(:,4) == 0);

S.File = File;
S.Candidates = [Data(r,[1 2]) ones(length(r),1)];
S.Query = [];
S.Discrepancy = ones(length(r),1);

xListCandidates(S);

%xDisplayCandidates(File,

