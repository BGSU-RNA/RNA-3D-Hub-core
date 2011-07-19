% zMarkRedundantChains(File) identifies redundant chains in File and sets File.Redundant(i,j) = 1 if nucleotides with indices i and j, from different chains, are redundant when these chains are aligned and superimposed

function [File,LongestChain] = zMarkRedundantChains(File, Verbose)

if nargin < 2,
  Verbose = 1;
end

DiscCutoff = 0.4;

for f = 1:length(File),
 File(f).Redundant = sparse(zeros(length(File(f).NT)));
 
 if length(File(f).NT) > 0,
  Chain = cat(2,File(f).NT.Chain);              % all chain identifiers
  U = unique(Chain);                            % unique chain identifiers

  if length(U) > 1,                             % more than one chain

    for u = 1:length(U),                        % loop through chains
      i = find(Chain == U(u));                  % indices of this chain
      bases{u} = cat(2,File(f).NT(i).Base);     % store bases of this chain
      indic{u} = i;                             % store indices of this chain
      leng(u)  = length(i);                     % store length of this chain
    end

    if Verbose > 0,
      fprintf('%s chains:\n', File(f).Filename);
      for u = 1:length(U),
        fprintf('Chain %s:  %s\n', U(u), bases{u});
      end
      fprintf('\n');
    end

    [y,i] = max(leng);
    LongestChain{f} = U(i);                     % chain with greatest length

    % ------------------------------------------- compare chains

    pro = zeros(length(U));
    d = ones(length(U));                        % default discrepancy btw chain

    for u = 1:length(U),                        
      for v = (u+1):length(U),

        % --------------------------------------- align sequences
        if max(leng(u),leng(v)) / min(leng(u),leng(v)) > 1.5,
          pro(u,v) = 0;                         % very different lengths
        else
          [matches,a,b,ss,tt] = dNeedlemanWunsch(bases{u},bases{v});
          e = find(bases{u}(a) == bases{v}(b)); % locations of agreement
          matches = length(e);                  % number of agreements
          a = a(e);
          b = b(e);
          pro(u,v) = matches/min(leng(u),leng(v));  % percent identity
        end

        pro(v,u) = pro(u,v);                    % symmetrize

%       if ((length(bases{u}) - matches < 4) || (pro > p)),

        % --------------------------------------- superimpose geometrically
        if pro(u,v) > 0,
          [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a),File(f),indic{v}(b));

          if Disc > DiscCutoff,                 % not such good superposition
            c = find(abs(A) <= 0.8);              % reasonably similar bases
            [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a(c)),File(f),indic{v}(b(c)));
          end

          d(u,v) = Disc;
          d(v,u) = d(u,v);
        end

      end
    end

d

    redundant = (pro >= 0.95) .* (d < 0.4);           % redundant chains
    redundant = redundant + eye(size(redundant));     % reflexivity
    redundant = redundant^10;                         % extend by transitivity

    if Verbose > 0,
      fprintf('Redundant chains for %4s: ', File(f).Filename);
    end

    for u = 1:length(U),
      for v = (u+1):length(U),
        if redundant(u,v) > 0,
          [matches,a,b,ss,tt] = dNeedlemanWunsch(bases{u},bases{v});
          e = find(bases{u}(a) == bases{v}(b));      % locations of agreement
          matches = length(e);
          a = a(e);
          b = b(e);

          File.Redundant(indic{u}(a),indic{v}(b)) = eye(matches);
          File.Redundant(indic{v}(b),indic{u}(a)) = eye(matches);

          if Verbose > 0,
            fprintf('%c%c ', U(u), U(v));
          end
        end
      end
    end

    if Verbose > 0,
      fprintf('\n');
    end

    if Verbose > 1,
      figure(1)
      clf
      zCircularDiagram(File(f),1);
      figure(2)
      clf
      spy(File(f).Redundant)
      drawnow
      pause
    end

  else
    LongestChain{f} = U(1);                     % chain with greatest length
  end
 end
end

