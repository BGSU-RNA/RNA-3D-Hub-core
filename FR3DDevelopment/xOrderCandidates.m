% xOrderCandidates orders candidates by similarity and displays the matrix of discrepancies between them

function [Search] = xOrderCandidates(File,Search,Level,UsingFull)

N      = Search.Query.NumNT;
Search = xMutualDiscrepancy(File,Search);      % compute discrepancy matrix
Done   = find(Search.DiscComputed);            % ones already computed

for i=1:length(Done),
  f = Search.Candidates(Done(i),N+1);          % file number
  b = '';
  for j = 1:min(4,N),
    b = [b File(f).NT(Search.Candidates(Done(i),j)).Base];
  end
  n = File(f).NT(Search.Candidates(Done(i),1)).Number;
  n = sprintf('%4s',n);
  if Search.Query.Geometric > 0,
      if isfield(Search,'AvgDisc'),
        d = sprintf('%6.4f',Search.AvgDisc(Done(i)));
      else
        d = sprintf('%6.4f',Search.Discrepancy((Done(i))));
      end
    else
      d = sprintf('%5d',Search.Discrepancy(Done(i))); % orig candidate number
    end
  Lab{i} = [b n ' ' File(f).Filename];
end

D = Search.Disc(Done,Done);                    % mutual distances to consider

[s,t] = size(D);

L = length(Search.Candidates(:,1));

ClusterNum = L+1;
CurrentCluster = 1:L;

for i = 1:L,
  Search.GroupLabel{i} = '          ';
end

figure(99)
clf
p = zOrderbySimilarity(D);
zGraphDistanceMatrix(D(p,p),Lab(p),15,0);
title('Table of discrepancies between candidates');

colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 0.8]);
colorbar('location','eastoutside');

%set(gcf,'Renderer','painters');
%saveas(gcf,[ 'Temp2.pdf'],'pdf');

n = length(p):-1:1;
r = Done(p(n));
q = Done(p);                                      % reverse order

% ----------------------------------- List and display results

fprintf('Candidates are ordered by similarity:\n');
fprintf('\n');

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(q,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(q);
S.Disc         = Search.Disc(q,q);
S.DiscComputed = Search.Disc(1,q);
S.File         = Search.File;
S.GroupLabel   = Search.GroupLabel(q);
if isfield(Search,'AvgDisc'),
  S.AvgDisc = Search.AvgDisc(q);
end

S.CandidateFilenames = Search.CandidateFilenames;

xListCandidates(S,Inf);                             % show on screen

xDisplayCandidates(File,S,Level+1,UsingFull);       % display, level 1
