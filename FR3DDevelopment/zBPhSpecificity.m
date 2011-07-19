% zBPhSpecificity reads Jesse's dataset on E.coli with number of interactions and potential base-protein interactions, and makes simple counts of base frequency when certain BPh interactions are observed

[n,t] = xlsread('Bases_from_Ec_and_Interactions_10_21-1.xls');

t = t(3:end,:);                               % remove first two rows

Filenames = {'2AVY','2AW4'};

 File = zAddNTData(Filenames);
% File = zAttachAlignment(File,1);
% xAnnotate

Data = [];
Spec = [];

for f = 1:length(File),
  BP(f).BP = File(f).BasePhosphate;
  for nn = 1:File(f).NumNT,
    BP(f).BP(nn,nn) = 0;                      % no self interactions
  end
end

c = 1;                                        % BPh counter

for r = 1:length(n(:,1)),                     % go through alignment
  f = find(ismember(lower(Filenames),lower(t{r,1})));      % get file number
  % fprintf('%d %s  %s\n', f, Filenames{f}', t{r,1})
  Data(r,1) = f;                                           % record file num
  switch f,
    case 1, i = zIndexLookup(File(f),num2str(n(r,1)),'A'); % look up indices
    case 2, i = zIndexLookup(File(f),num2str(n(r,1)),'B');
  end
  i = i(1);
  Data(r,2) = i;                                % index of this nucleotide
  Data(r,17) = File(f).NT(i).Code;              % which base this is

  Data(r,8) = sum((abs(File(f).Edge(i,:)) > 20) .* (abs(File(f).Edge(i,:)) < 24));                                              % does it stack?

  Data(r,10) = sum(((BP(f).BP(i,:)) > 0) .* ((BP(f).BP(i,:)) < 30));
  Data(r,12) = sum(((BP(f).BP(:,i)) > 0) .* ((BP(f).BP(:,i)) < 30));

  Data(r,15) = 0;
  Data(r,16) = 0;

  Data(r,13) = sum((fix(abs(File(f).Edge(i,:))) ==1));
  Data(r,14) = sum((abs(File(f).Edge(i,:)) >= 2) .* (abs(File(f).Edge(i,:)) < 15));

  counts = n(r,4:7);
  percs  = 0.00001 + counts / sum(counts);        % number of A, C, G, U

  Data(r,11) = - percs * log2(percs)';            % entropy of the distribution
  Data(r,3) = 100*counts(File(f).NT(i).Code)/sum(counts); % conservation

  File(f).Conservation(i) = Data(r,3);

  bph = find(((BP(f).BP(i,:)) > 0) .* ((BP(f).BP(i,:)) < 30));

  for b = bph,                                  % run through BPh inter
    Spec(c,1) = BP(f).BP(i,b);                  % record BPh type
    Spec(c,2) = Data(r,13);                     % # cWW made
    Spec(c,3) = Data(r,14);                     % # non-cWW made
    Spec(c,4:7) = counts;                       % counts from sequences
    Spec(c,8)   = File(f).NT(i).Code;           % base making the BPh
    c = c + 1;
  end

  if isnan(Data(r,11)),
    counts
    Data(r,3)
  end
end

N = length(Spec(:,1));                          % number of instances
for c = 1:N,
  Spec(c,9:12) = Spec(c,4:7) / sum(Spec(c,4:7)); % normalize
  Spec(c,3)    = min(1,Spec(c,3));               % 0 if no non-cWW pair
  Spec(c,13)   = min(1,Spec(c,2) + Spec(c,3));   % 0 if no pair, 1 if pair
end

Spec = sortrows(Spec);

for c = 1:4,
  j = find(Spec(:,8) == c);                     % select by base
  ByBase(c,1:4) = mean(Spec(j,9:12));
  ByBase(c,5) = length(j);
end

ByBase

for bp = 0:1,
  for c = 1:4,
    j = find((Spec(:,8) == c) .* (Spec(:,2) == bp));   % select by base
    BycWWPairByBase(c,1:4) = mean(Spec(j,9:12));
    BycWWPairByBase(c,5) = length(j);
  end
  BycWWPairByBase
end

for bp = 0:1,
  for c = 1:4,
    j = find((Spec(:,8) == c) .* (Spec(:,3) == bp));   % select by base
    BynoncWWPairByBase(c,1:4) = mean(Spec(j,9:12));
    BynoncWWPairByBase(c,5) = length(j);
  end
  BynoncWWPairByBase
end

for bp = 0:1,
  for c = 1:4,
    j = find((Spec(:,8) == c) .* (Spec(:,13) == bp));   % select by base
    ByPairByBase(c,1:4) = mean(Spec(j,9:12));
    ByPairByBase(c,5) = length(j);
  end
  ByPairByBase
end




% ----------------------------------------------- select by active edge
ae    = {'SW',' W',' H',' H',' W',' H',' H',' H',' H',' S',' W',' W',' W',' H',' W',' H',' H',' H',' W'};

edge  = [4 1 2 2 1 2 2 2 2 3 1 1 1 2 1 2 2 2 1]; % 1 = W, 2 = H, 3 = S, 4 = SW
aelet = {'W','H','S','T'};
aelet(edge)
Letter = 'ACGU';

for c = 1:4,
  for e = 1:4,
    j = find((edge(Spec(:,1))' == e) .* (Spec(:,8) == c)); % select edge, base
    ByBaseByEdge(e,1:4) = mean(Spec(j,9:12));
    ByBaseByEdge(e,5) = length(j);
  end

  fprintf('Base %s:\n', Letter(c));
  ByBaseByEdge

end


for bph = 1:19,
  j = find(Spec(:,1) == bph);                     % select by base
  ByBPh(bph,1:4) = mean(Spec(j,9:12));
  ByBPh(bph,5)   = length(j);
end

ByBPh

