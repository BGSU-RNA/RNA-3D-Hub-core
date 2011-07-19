% zUniqueChains re-orders the nucleotides in File to put the largest chain first

% File = zAddNTData({'2BU1','409D','2AWE','1N7A','1QCU'});
% File = zStoreO3(File);

function [File,Indices] = zOrderChains(File, Verbose)

if nargin < 2,
  Verbose = 0;
end

p = 0.95;

%figure(1)
%clf

for f = 1:length(File),

  clear pro
  clear bases
  clear indic
  clear d
  inter = [];

  Chain = cat(2,File(f).NT.Chain);

  U = unique(Chain);
  clear bases

 if length(U) > 1,

  for u = 1:length(U),
    i = find(Chain == U(u));
    bases{u} = cat(2,File(f).NT(i).Base);
    indic{u} = i;
  end

%  bases

  for u = 1:length(U),
    for v = (u+1):length(U),
      [matches,a,b,ss,tt] = zNeedlemanWunsch(bases{u},bases{v});
      pro(u,v) = matches/min(length(bases{u}),length(bases{v})); % percent agreement

      lu = length(bases{u});
      lv = length(bases{v});

      if max(lu,lv) / min(lu,lv) > 1.5,
        pro(u,v) = 0;
      end

      pro(v,u) = pro(u,v);
%      if ((length(bases{u}) - matches < 4) || (pro > p)),
      if pro(u,v) > 0,
        d(u,v) = xDiscrepancy(File(f),indic{u}(a),File(f),indic{v}(b));
        d(v,u) = d(u,v);
      end

      e = abs(File(f).Edge(indic{u},indic{v}));
      inter(u,v) = sum(sum((e>0).*(e<30)));
      inter(v,u) = inter(u,v);

%      zShowInteractionTable(File(f),[indic{u},indic{v}]);

    end
  end

  inter = full(inter);

%  pro
%  d
%  inter

%  figure(1)
%  subplot(3,2,f)
%  spy(File(f).Edge)

  % distinct is a relationship between two chains

  keeptogether = full((inter > 0) + (pro <= 0.95));

  k = keeptogether^20;

  if Verbose > 0,
    fprintf('zUniqueChains: File %4s has chains %10s.  Keeping ', File(f).Filename, U);
  end

  Indices{f} = [];
  for i = 1:length(U),
    if k(1,i) > 0,
      Indices{f} = [Indices{f} indic{i}];
      if Verbose > 0,
        fprintf('%c', U(i));
      end
    end
  end

  if Verbose > 0,
    fprintf('\n');
  end


%  clf
%  VP.Sugar = 1;
%  zDisplayNT(File(f),Indices{f},VP);

%  pause
 else
   Indices{f} = 1:length(File(f).NT);
 end
end

