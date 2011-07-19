
% File = zAddNTData({'2avy','1j5e'});
Verbose = 1;
a = 0;

% -------------------------------------- Load JesseRyan composite alignment

load JesseRyan16SComposite

a = a+1;
Al(a).ModelStructure = JRComposite2AVY;
Al(a).InferStructure = JRComposite1J5E;
Al(a).Name           = 'Jesse-Ryan composite';

% -------------------------------------- Load instances from a search

load 2009-07-21_14_13_55-2avy_triples

Cand = Search.Candidates;
[L,N] = size(Cand);
N = N - 1;

% -------------------------------------- Loop through candidates, list
  fprintf('          ');
  for n = 1:N,
    fprintf('%1s %5d ', ' ', n);
  end
  for n = 1:N,
    for m = (n+1):N,
      fprintf('%6s ', [num2str(n) '-' num2str(m)]);
    end
  end
  for n = 1:N,
    fprintf('%1s %5d ', ' ', n);
  end
  for n = 1:N,
    for m = (n+1):N,
      fprintf('%6s ', [num2str(n) '-' num2str(m)]);
    end
  end
  fprintf('\n');

clear j

for c = 1:L,
  fprintf('Cand %3d: ', c);
  for n = 1:N,
    NT = File(1).NT(Cand(c,n));
    fprintf('%1s %5s ', NT.Base, NT.Number);
  end
  for n = 1:N,
    for m = (n+1):N,
      fprintf('%6s ', zEdgeText(File(1).Edge(Cand(c,n),Cand(c,m)),1));
    end
  end
  for n = 1:N,
    k = find(Al(a).ModelStructure == Cand(c,n));
    j{n} = Al(a).InferStructure(k);
    if ~isempty(j{n}),
      NT = File(2).NT(j{n}(1));
      fprintf('%1s %5s ', NT.Base, NT.Number);
    else
      fprintf('%1s %5s ', ' ', '     ');
    end
  end
  for n = 1:N,
    for m = (n+1):N,
      if ~isempty(j{n}) && ~ isempty(j{m}),
        fprintf('%6s ', zEdgeText(File(2).Edge(j{n},j{m}),1));
      end
    end
  end
  fprintf('\n');
end

