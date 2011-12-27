% xPairsAndIntersections.m loads a PDB file, selects a query motif randomly, determines FR3D search parameters, does pairwise screening, and intersects the pairwise results

Motif   = 2;
Verbose = 1;
clear Query

switch Motif
case 1,
  FN = '1J1U';
  Query.Description    = 'FR3D test query';
  Query.Filename       = FN;
  m = 4;                        % number of nucleotides in query
  Query.IndicesManual  = 10:(10+m-1);                   % indices
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {FN};
  Query.Name           = 'Test';
  Query.NTList         = {1,1,1,1};
case 2,
  FN = '1S72';
  m = 5;
  Query.Description    = '1S72 Sarcin core';
  Query.Filename       = FN;
  Query.IndicesManual  = [1251 1252 1253 898 899];       % indices
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {FN};
  Query.Name           = 'Sarcin core';
  Query.Geometric      = 1;
  Query.Mask           = 'NNNNN';
  Query.NTList         = {1,1,1,1,1};            % fake list, correct length
end

Filenames = {FN};
File = zAddNTData(FN);

UsingLibrary = 1;

Query = xConstructQuery(Query,File);
Query.DistCutoff = 30;

c = cat(1,File(1).NT(1:File(1).NumNT).Center);
File(1).Distance = zMutualDistance(c,Query.DistCutoff); 

Codes = cat(1,File.NT(:).Code); 

Model = Query;

% -------------------------------- Create file name from PDB ID and NTs
OF = ['Database' filesep FN];
for i = 1:length(Query.NT),
  OF = [OF '_' Query.NT(i).Number];
end

% -------------------------------- Write canditates
OFN = [OF '_Pairs.txt'];
fid = fopen(OFN,'w');

for i=2:Model.NumNT                % loop through model nucleotides
  for j=1:(i-1)
    PS{i,j} = xPairwiseScreen(File,Codes,Model,i,j);
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

% -------------------------------- Write canditates
OFN = [OF '_Candidates.txt'];
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

