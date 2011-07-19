% zPhosphateSubstitutions reads Jesse's alignments of 16S and base-phosphate interactions and substitution data from sequence alignments, then makes graphs for each base-phosphate interaction category of the pattern of substitutions for the PHOSPHATE ACCEPTOR

z = 2;

if z == 1,
  [n,t] = xlsread('16S_Ec_Tt_BPh_Aln_12_6_08_Seq_Data.xls');
  t = t(3:end,:);
  [nn,tt] = xlsread('23S_Ec_Tt_BPh_Aln_12_6_08_Seq_Data.xls');
  tt = tt(3:end,:);
else
  [n,t] = xlsread('16S_Ec_Tt_BPh_Aln_12_23_08.xls');
  t = t(3:end,:);
  [nn,tt] = xlsread('23S_Ec_Tt_BPh_Aln_12_23_08.xls');
  tt = tt(3:end,:);
  [a,b] = size(n);
  [aa,bb] = size(t);
  if a ~= aa,
    fprintf('Fix the first file; n and t have different numbers of rows!');
  end
end

tt = [t; tt];
nn = [n; nn];

[a,b] = size(nn);

for i = 1:a,
  for j = 1:b,
    if isnan(nn(i,j)),
      nn(i,j) = 0;
    end
  end
end

tt(:,10) = strrep(tt(:,10),'8bBPh','7BPh');
tt(:,10) = strrep(tt(:,10),'4bBPh','4BPh');

if 0 > 1,
  FN = unique(tt(:,4));
  for f = 1:length(FN),
    FNL(f) = length(FN{f,1});
  end
  FN = FN(find(FNL));
  File = zAddNTData(FN);
else
  FN = {'2AVY','2AW4','2AW7','2AWB'};
  File = zAddNTData(FN);
  Chain = {'A','B','A','B'};
end

for y = 1:length(tt(:,1)),
  f = find(ismember(FN,tt{y,4}));
  if ~isempty(f),
    [i,ac] = zIndexLookup(File(f),num2str(nn(y,6)),Chain{f});
    e = abs(fix(File(f).Edge(i(1),:)));
    if sum( (e >= 2) .* (e < 15) ) > 0,                % makes a non-cWW
      if sum(e == 1) > 0,     % makes a cWW
        inter(y) = 1;
      else
        inter(y) = 2;      
      end
    else
      if sum(e == 1) > 0,     % makes a cWW
        inter(y) = 3;
      else
        inter(y) = 4;      
      end
    end
  end
end


A = [0 1];
C = [1 1];
U = [1 0];
G = [0 0];

Letter = {'A','C','G','U'};
color = [[1 0 0]; [0.9 0.9 0]; [0 0.9 0]; [0 0 1]];
types = {'cWW and non-cWW','non-cWW only','cWW only','no basepair'};
figure(z)
clf
for t = 1:4,

  subplot(2,2,t);

  % ----------------------- aligned positions plotted as dots
  p = find(inter == t);                % rows with OK BPh type

  i = zeros(length(p),1);
  for j = 1:length(p),
    if ~isempty(tt{p(j),19}),
      i(j) = 1;
    end
  end
  p = p(find(i));

  % Use counts for the base phosphate acceptor
  counts = sum(nn(p,22:25),2);                          % sequence counts
  L = nn(p,29)*A + nn(p,30)*C + nn(p,31)*G + nn(p,32)*U ;  % sequence counts

  L = L ./ (counts * [1 1]);                            % weighted average
  for i = 1:4,
    j = find(cat(1,tt{p,7}) == Letter{i});              % base in 3D structure
    scatter(L(j,1),L(j,2),10,color(i,:),'filled');      % color differently
    hold on
  end
  
  pp = length(p);

  title([num2str(pp) ' instances of ' types{t}],'fontsize',8);
  axis([0 1 0 1]);
  axis off
  s = 0.07;
  text(-s,1,'A','HorizontalAlignment','center','Color',color(1,:));
  text(1+s,1,'C','HorizontalAlignment','center','Color',color(2,:));
  text(1+s,0,'U','HorizontalAlignment','center','Color',color(4,:));
  text(-s,0,'G','HorizontalAlignment','center','Color',color(3,:));
  plot([0 1 1 0 0], [0 0 1 1 0], 'k');
  axis square
end

saveas(gcf,['Phosphate Interactions\PhosphateSubstitutions_phosphate.pdf'],'pdf');
