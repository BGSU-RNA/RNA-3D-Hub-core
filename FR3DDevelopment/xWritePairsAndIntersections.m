% xPairsAndIntersections.m loads a PDB file, selects a query motif randomly, determines FR3D search parameters, does pairwise screening, and intersects the pairwise results

FN = '1J1U';
Verbose = 1;
m = 4;                        % number of nucleotides in query

clear Query

% ---------------------------------------------------------------------------

File = zAddNTData(FN);
i = 1;

Filenames = {FN};

  Query.Description    = 'FR3D test query';
  Query.Filename       = FN;
  Query.IndicesManual  = 10:(10+m-1);                   % indices
% Query.NTList ?????
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {FN};
  Query.Name           = 'Test';

UsingLibrary = 1;

Query = xConstructQuery(Query,File);

Query.DistCutoff = 30;

    c = cat(1,File(i).NT(1:File(i).NumNT).Center);
    File(i).Distance = zMutualDistance(c,Query.DistCutoff); 

Codes = cat(1,File.NT(:).Code); 

Model = Query;
f = 1;

OFN = ['Database' filesep FN];
for i = 1:length(Query.NT),
  OFN = [OFN '_' Query.NT(i).Number];
end
OFN = [OFN '_Pairs.txt'];

fid = fopen(OFN,'w');

for i=2:Model.NumNT                % loop through model nucleotides
  for j=1:(i-1)
    PS{i,j} = xPairwiseScreen(File(f),Codes,Model,i,j);
    PS{j,i} = PS{i,j}';
    NNZ(i,j) = nnz(PS{i,j});              % number of non-zero entries
    NNZ(j,i) = NNZ(i,j);

    fprintf(fid,'%d,%d',j,i);
    [a,b] = find(PS{i,j});
    for c = 1:length(a),
      fprintf(fid,',%s,%s',File.NT(b(c)).Number,File.NT(a(c)).Number);
    end
    fprintf(fid,'\n');
  end
end
fclose(fid);

Query.SSCutoff = Inf*Query.SSCutoff;

Candidates = xFindCandidates(File,Query,Verbose);  % screen for candidates

OFN = ['Database' filesep FN];
for i = 1:length(Query.NT),
  OFN = [OFN '_' Query.NT(i).Number];
end
OFN = [OFN '_Candidates.txt'];
fid = fopen(OFN,'w');

for c = 1:length(Candidates(:,1)),
  for n = 1:length(Query.NT),
    if n > 1,
      fprintf(fid,',%s',File.NT(Candidates(c,n)).Number);
    else
      fprintf(fid,'%s',File.NT(Candidates(c,n)).Number);
    end
  end
  fprintf(fid,'\n');
end
fclose(fid);

break

if (Query.Geometric > 0),

  % --------- Compute square of distance difference from model

  d = (d - Query.Distance(p,q)).^2;   % squared difference in distances

  d = d + 0.00000001 * (d == 0);       % avoid rejecting model; make d nonzero

  % --------- Impose upper limit on distance differences; 2-nucleotide cutoff


  if Query.NumNT > 2,
    Wp = Query.LocWeight(p);
    Wq = Query.LocWeight(q);
    MaxD = (Wp + Wq) * (Query.NumNT * Query.DiscCutoff)^2 / (Wp * Wq);
  else
    Wp = 1;
    Wq = 1;
    MaxD = (Query.NumNT * Query.DiscCutoff)^2;
  end

  if isfield(Query,'Flex'),
    MaxD = max(MaxD,Query.Flex(p,q)^2);  % allow larger distance if desired
  end

  k = find(d <= MaxD);            % keep ones with small difference from model

  i = i(k);
  j = j(k);
  d = d(k) * Wp * Wq;

end
