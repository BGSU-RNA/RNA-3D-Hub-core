% xPairSubstitutions(File,f,i,j) does a FR3D geometric search in File using the pair in file f, indices i and j, as the query.  It returns a list of the codes of what is found along with the discrepancy (in D1) and IDI (in D) from the query pair

% If i and j are nucleotides, then f is the code for the interaction between them

%Filenames = {'2avy','2aw4'};
%File = zAddNTData(Filenames);

%f = 1;

%i = zIndexLookup(File(f),'55');
%j = zIndexLookup(File(f),'56');

function [D,D1] = xPairSubstitutions(File,f,i,j)

if ~strcmp('double',class(i)),
  F.Filename = 'synthetic';
  F.NT(1) = i;
  F.NT(2) = j;
  F.Edge(1,2) = f;
  F.Edge(2,1) = -f;
  f = 1;
  i = 1;
  j = 2;
else
  F = File(f);
end

for a = 1:length(File),
  Filenames{a} = File(a).Filename;
end

Verbose = 0;

Query.Name           = 'Pairwise interaction search';
Query.Description    = 'Pairwise interaction search';
Query.Filename       = F.Filename;
Query.NTList         = {F.NT(i).Number F.NT(j).Number};
Query.ChainList      = {F.NT(i).Chain F.NT(j).Chain}; 
Query.Edges{1,2}     = [zEdgeText(F.Edge(i,j)) ' n' zEdgeText(F.Edge(i,j))];
Query.DiscCutoff     = 0.9;      

Query = xConstructQuery(Query,F);

UsingLibrary = 1;                  

xFR3DSearch

L = length(Candidates(:,1));        % number of candidates found

D = zeros(L,3);

for c = 1:L,
  D(c,1) = File(Candidates(c,3)).NT(Candidates(c,1)).Code;
  D(c,2) = File(Candidates(c,3)).NT(Candidates(c,2)).Code;
  D(c,3) = Discrepancy(c);
end

D1 = D;

%Search
%Search.Query

[Discrepancy, Candidates] = xRankCandidatesIDI(File,Search.Query,Candidates,Verbose);

D = zeros(L,3);

for c = 1:L,
  D(c,1) = File(Candidates(c,3)).NT(Candidates(c,1)).Code;
  D(c,2) = File(Candidates(c,3)).NT(Candidates(c,2)).Code;
  D(c,3) = Discrepancy(c);
end

% xListCandidates(Search,20);
