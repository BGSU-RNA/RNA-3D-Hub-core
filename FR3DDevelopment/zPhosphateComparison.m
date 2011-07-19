% zPhosphateComparison takes phosphate classifications from multiple files and calculates average phosphate interaction values

% cat(2,File(3).NT(:).Chain)

Set = 8;

switch Set,

case 1,
Representative = '1s72';
% File = zAddNTData('1s72_equiv');             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison1s72 D
% load PhosphateComparison1s72
OK = '0A';                                    % use these chains with 1s72

case 2,
Representative = '2avy';
% File = zAddNTData({'2avy','2aw7'});
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison2avy D
% load PhosphateComparison2avy
OK = 'A';

case 3,
Representative = '1s72Self';
% File = zAddNTData('1s72_equiv');             % load data
% [File,D] = zPhosphateInteractions(File,2,1); % analyze interactions
% save PhosphateComparison1s72Self D
% load PhosphateComparison1s72Self
OK = '0A';                                    % use these chains with 1s72

case 4,
Representative = '2j01';
% File = zAddNTData({'2j01','2j03','1vsa'});             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison2j01 D
% load PhosphateComparison2j01
OK = 'Aw';                                    % use these chains with 1s72

case 6,
Representative = '1j5e';
% File = zAddNTData({'1j5e','2j00','2j02'});             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison2j01 D
% load PhosphateComparison2j01
OK = 'Aw';                                    % use these chains with 1s72

case 7,
Representative = '2j01-2aw4';
% File = zAddNTData({'2j01','2aw4'});             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison2j01-2aw4 D
load PhosphateComparison2j01-2aw4
OK = 'BA';                                    % use these chains
OK = 'AB';                                    % use these chains

case 8,
Representative = '2aw4-2awb';
% File = zAddNTData({'2aw4','2awb'});             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparison2aw4-2awb D
load PhosphateComparison2aw4-2awb
OK = 'BB';                                    % use these chains

end

% --------------------------------------------- preliminary data work

i = find(D(:,17) == 1);                        % use best oxygen only
D = D(i,:);

D(:,24) = mod(D(:,5),100);                      % don't distinguish near BPh

D = sortrows(D,5);                              % sort by classification

%D = sortrows(D,[21 22 1]);

% --------------------------------------------- restrict to listed chains

i = zeros(size(D(:,1)));

for m = 1:length(D(:,1)),
  if (File(D(m,1)).NT(D(m,2)).Chain == OK(D(m,1))) && (File(D(m,1)).NT(D(m,3)).Chain == OK(D(m,1))),    % specify which structure uses which chain
%  if any(File(D(m,1)).NT(D(m,2)).Chain == OK) && any(File(D(m,1)).NT(D(m,3)).Chain == OK),
    i(m) = 1;
  end
end

D = D(find(i),:);

% --------------------------------------------- Nucleotide-nucleotide inter

[s,t] = size(D);

n1n2 = unique(D(:,[21 22]), 'rows');        % unique pairs of nucleotides

L = zeros(length(n1n2),119);

cc = 1;

for m = 1:length(n1n2(:,1)),                   % loop through unique inters
%  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)) .* (D(:,23) == n1n2(m,3)));
  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)));

  [a,b] = sort(D(i,1));                         % sort by file number
  i = i(b);

  if min(D(i,4)) < max(D(i,4)),
%     [n1n2(m,:) min(D(i,4)) max(D(i,4))]
%     D(i,:)                                     % check if bases are the same
  end

  PrevFile = 0;

flag = 0;

  for k = 1:length(i),
    if D(i(k),1) ~= PrevFile,
      L(m,D(i(k),23)) = L(m,D(i(k),23)) + 1;

if D(i(k),23) < 100,
fprintf('File %4s Base %4d Phosphate %4d Interaction %5s Number %3d\n', File(D(i(k),1)).Filename, D(i(k),21), D(i(k),22), zBasePhosphateText(D(i(k),23)),cc);
flag = 1;
end

      PrevFile = D(i(k),1);
    end
  end

if flag == 1,
  cc = cc + 1;
end


end

L = sortrows(L,1:119);

A = [0 1];
C = [1 1];
G = [1 0];
U = [0 0];

figure(1)
clf
  Letter = {'3BPh','4BPh','4BPh-bif','5BPh'};
  M = L(:,[11 12 19 13]);
  i = find(sum(M,2));                  % rows having a true BPh of these types
  M = M(i,:);
  nn = M;
  counts = sum(M,2);
  LL = nn(:,1)*A + nn(:,2)*C + nn(:,3)*G + nn(:,4)*U ;
  LL = LL ./ (counts * [1 1]);
  scatter(LL(:,1),LL(:,2),10,counts,'filled');
  hold on
  axis([0 1 0 1]);
  axis off
  colormap('default');
  map = colormap;
  colormap(map(1:58,:));
  colorbar('southoutside');
  s = 0.07;
  text(-s,1,Letter{1},'HorizontalAlignment','center');
  text(1+s,1,Letter{2},'HorizontalAlignment','center');
  text(1+s,0,Letter{3},'HorizontalAlignment','center');
  text(-s,0,Letter{4},'HorizontalAlignment','center');
  plot([0 1 1 0 0], [0 0 1 1 0], 'k');
  title(['Variation between ' Letter{1} ' ' Letter{2} ' ' Letter{3} ' ' Letter{4} ' in ' num2str(length(File)) ' ' Representative ' equivalents']);

  saveas(gcf,['Phosphate Interactions\PhosphateChanges' Representative '_G_WC_edge.pdf'],'pdf');

figure(2)
clf
subplot(2,1,1)
hist(nonzeros(sum(M,2)),1:length(File));
nbin = hist(nonzeros(sum(M,2)),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar G BPh interaction');

subplot(2,1,2)
hist(nonzeros(sum(M,2)+sum(L(i,[111 112 119 113]),2)),1:length(File));
nbin = hist(nonzeros(sum(M,2)+sum(L(i,[111 112 119 113]),2)),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar G BPh or nBPh interaction');

  saveas(gcf,['Phosphate Interactions\PhosphateChanges' Representative '_G_WC_edge_Hist.pdf'],'pdf');

figure(3)
clf
  Letter = {'7BPh','8BPh','8BPh-bif','9BPh'};
  M = L(:,[6 7 18 8]);
  i = find(sum(M,2));                  % rows having a true BPh of these types
  M = M(i,:);
  nn = M;
  counts = sum(M,2);
  LL = nn(:,1)*A + nn(:,2)*C + nn(:,3)*G + nn(:,4)*U ;
  LL = LL ./ (counts * [1 1]);
  scatter(LL(:,1),LL(:,2),10,counts,'filled');
  hold on
  axis([0 1 0 1]);
  axis off
  colormap('default');
  map = colormap;
  colormap(map(1:58,:));
  colorbar('southoutside');
  s = 0.07;
  text(-s,1,Letter{1},'HorizontalAlignment','center');
  text(1+s,1,Letter{2},'HorizontalAlignment','center');
  text(1+s,0,Letter{3},'HorizontalAlignment','center');
  text(-s,0,Letter{4},'HorizontalAlignment','center');
  plot([0 1 1 0 0], [0 0 1 1 0], 'k');
  title(['Variation between ' Letter{1} ' ' Letter{2} ' ' Letter{3} ' ' Letter{4} ' in ' num2str(length(File)) ' ' Representative ' equivalents']);

  saveas(gcf,['Phosphate Interactions\PhosphateChanges' Representative '_C_Hoogsteen_edge.pdf'],'pdf');

figure(4)
clf
subplot(2,1,1)
hist(nonzeros(sum(M,2)),1:length(File));
nbin = hist(nonzeros(sum(M,2)),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar C BPh interaction');

subplot(2,1,2)
hist(nonzeros(sum(M,2)+sum(L(i,[106 107 118 108]),2)),1:length(File));
nbin = hist(nonzeros(sum(M,2)+sum(L(i,[106 107 118 108]),2)),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar C BPh or nBPh interaction');

saveas(gcf,['Phosphate Interactions\PhosphateChanges' Representative '_C_Hoogsteen_edge_Hist.pdf'],'pdf');

% --------------------------------------------- Histograms over all phosphates

[s,t] = size(L);

N = L;
N(:,11) = sum(L(:,[11 12 19 13]),2);          % roll together G WC inters
N(:,[12 19 13]) = zeros(s,3);
N(:,111) = sum(L(:,[111 112 119 113]),2);          % roll together G WC inters
N(:,[112 119 113]) = zeros(s,3);
N(:, 6) = sum(L(:,[6 7 18 8]),2);             % roll together C WC inters
N(:,[7 18 8]) = zeros(s,3);
N(:,106) = sum(L(:,[106 107 118 108]),2);             % roll together C WC inters
N(:,[107 118 108]) = zeros(s,3);

clear Y

c = 1;
for r = 1:s,
  p = find(N(r,1:19));                        % where are the nonzeros
  if length(p) > 1,
    p
  elseif length(p) == 1,                      % found a base-phosphate
    Y(c,1) = N(r,p(1));
    Y(c,2) = N(r,p(1)) + N(r,p(1)+100);        % this class and near    
    c = c + 1;
  end
end

fprintf('Found %d distinct base-phosphate interactions in %s and its equivalents\n', c-1, Representative);

figure(5)
clf

subplot(2,1,1)
nbin = hist(Y(:,1),1:length(File))
hist(Y(:,1),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar BPh interaction');

subplot(2,1,2)
nbin = hist(Y(:,2),1:length(File))
hist(Y(:,2),1:length(File));
axis([0.5 length(File)+0.5 0 1.1*max(nbin)]);
title('Number of structures having a similar BPh or nBPh interaction');

saveas(gcf,['Phosphate Interactions\PhosphateChanges' Representative '_All.pdf'],'pdf');

break

% --------------------------------------------- Individual interactions

[s,t] = size(D);


n1n2 = unique(D(:,[21 22 24]), 'rows');        % unique BPh interactions

clear E

L = zeros(length(n1n2),119);

for m = 1:length(n1n2(:,1)),                   % loop through unique inters
  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)) .* (D(:,24) == n1n2(m,3)));

%  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)) );

%length(i)
%D(i,[2 3 5])
%pause  

  if min(D(i,4)) < max(D(i,4)),
    D(i,:)                                     % check if bases are the same
  end

  if length(i) > 1,
    E(m,1:t) = mean(D(i,:));                   % average parameters
  else
    E(m,1:t) = D(i,:);                         % unless there is just one
  end
  E(m,1)   = D(i(1),1);                        % file number of first
  E(m,2)   = D(i(1),2);
  E(m,3)   = D(i(1),3);
  E(m,5)   = D(i(1),5);

  E(m,t+1) = length(i);                        % store number of such inters
  E(m,t+2) = std(D(i,9));                      % variation in bond length
  E(m,t+3) = std(D(i,8));                      % variation in angle

  for k = 1:length(i),
    L(m,D(i(k),5)) = L(m,D(i(k),5)) + 1;
  end
end

L = sortrows(L,1:119);

% --------------------------------------- remove cases of only near BPh

i = find(E(:,5) < 100);                          % first happens to be true BP
E = E(i,:);                                      % use these only


figure(1)
clf
n = hist(E(:,t+1),0.5:1:(length(File)-0.5));
hist(E(:,t+1),0.5:1:(length(File)-0.5));
xlabel('Number of structures with each base-phosphate interaction (near or true)');
axis([0 length(File) 0 max(n)*1.1]);

figure(2)
clf
scatter(E(:,9),E(:,8),6,E(:,t+1),'filled');
caxis([1 1.1*max(D(:,1))]);
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by number of structures having the same interaction');

figure(3)
clf
scatter(E(:,9),E(:,8),6,E(:,t+2),'filled');
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by standard deviation of bond length');

figure(4)
clf
scatter(E(:,9),E(:,8),6,E(:,t+3),'filled');
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by standard deviation of angle');

figure(5)
clf
plot(E(:,24)+0.5*rand(size(E(:,24)))-0.25,E(:,t+1),'.');
title('Which interactions are conserved between structures');
xlabel('Internal interaction code');
ylabel('Number of structures');

figure(6)
clf

n = hist(E(:,t+1),0.5:1:(length(File)-0.5));
hist(E(:,t+1),0.5:1:(length(File)-0.5));
xlabel('Number of structures with each base-phosphate interaction');
axis([0 length(File) 0 max(n)*1.1]);

break

% ------------------------------------- Display nucleotides in all structures
% ------------------------------------- when one has a true Base-Phosphate

% -------- Note: with 1s72_equiv, the large chain is 0 in some, A in others
% -------- The small chain is 9 or B.  So it is hard to find all the 
% -------- nucleotides in all structures.
% -------- Above, we include only chain 0 or A
% -------- Below, we look up the index using only nucleotide number

i = find(E(:,5) < 100);                          % first happens to be true BP
E = E(i,:);                                      % use these only

E = sortrows(E,[t+1 21 22]);                     % sort by number of inters

for r = 1:length(E(:,1)),                        % loop through instances
  Cand = [];

  for f = 1:length(File),                        % loop through files
    i = zIndexLookup(File(f),File(E(r,1)).NT(E(r,2)).Number,'0',0);
    i = [i zIndexLookup(File(f),File(E(r,1)).NT(E(r,2)).Number,'A',0)];

    j = zIndexLookup(File(f),File(E(r,1)).NT(E(r,3)).Number,'0',0);
    j = [j zIndexLookup(File(f),File(E(r,1)).NT(E(r,3)).Number,'A',0)];
    if ~isempty(i) && ~isempty(j),
      Cand = [Cand; [i j f]];
    end
  end

  xDisplayCandidates(File,Cand);
end


break

for f = 1:length(File),
  [b,i,j] = unique(cat(1,File(f).NT.Chain));
  fprintf('File %s has ', File(f).Filename);
  for m = 1:length(b),
    fprintf('%d in chain %s, ', sum(j==m), b(m));
  end
  fprintf('\n');
end

          % Columns of D
          % 1  file number
          % 2  index of base
          % 3  index of nucleotide using phosphate
          % 4  code of base
          % 5  classification number for this massive-oxygen pair
          % 6  which massive atom is interacting
          % 7  which oxygen is interacting
          % 8  angle of interaction, in degrees
          % 9  distance from massive atom to oxygen, in Angstroms
          %10  displacement of oxygen atom relative to C1' of base
          %13  displacement of phophorus atom relative to C1' of base
          %16  distance between centers of the two bases
          %17  1 if this is the best oxygen for this hydrogen, 0 otherwise
          %18  approximate quality of the hydrogen bond
          %19  angle made by massive, hydrogen, phosphorus
          %20  distance from massive to phosphorus
          %21  nucleotide number of base
          %22  nucleotide number of phosphate donor
          %23  classification of this interaction
          %24  interaction (without distinguishing true or near)