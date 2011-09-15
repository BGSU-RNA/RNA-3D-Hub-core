
function [Sig,RSig] = pUpdateSignature(MN)

load(['MotifLibrary' filesep MN '.mat']);

    [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
    N = N - 1;                              % number of nucleotides
    clear NewFile
    for ff = 1:length(Search.File),
      F = Search.File(ff);
      if ~isempty(F.NT),
        F.LooseCoplanar = sparse(F.NumNT,F.NumNT);
        NewFile(ff) = F;
      end
    end
    Search.File = NewFile;

CL = zClassLimits;                              % read ClassLimits matrix

if exist('PairExemplars.mat','file') > 0,
  load('PairExemplars','Exemplar');
else
  Exemplar = [];
end

    for c = 1:length(Search.Candidates(:,1)),
      ff = Search.Candidates(c,N+1);
      i = Search.Candidates(c,1:N);

      for a = 1:length(i),
        for b = (a+1):length(i),
          if Search.File(ff).Edge(i(a),i(b)) ~= 0,
            NT1 = Search.File(ff).NT(i(a));
            NT2 = Search.File(ff).NT(i(b));
            Pair = zLooseCoplanar(NT1,NT2,CL,Exemplar);
            Search.File(ff).LooseCoplanar(i(a),i(b)) = Pair.Coplanar;
            Search.File(ff).LooseCoplanar(i(b),i(a)) = Pair.Coplanar;
          end
        end 
      end
    end

[Edge,BPh,BR,Search] = pConsensusInteractions(Search); % consensus interactions

%full(Edge)

Sig = zMotifSignature(Edge,2,1,1);

index = find(diag(fix(abs(Edge)),1)==1);
Truncate = index+1;
    per = [Truncate:N 1:(Truncate-1)];
RSig = zMotifSignature(Edge(per,per),2,1,1);
