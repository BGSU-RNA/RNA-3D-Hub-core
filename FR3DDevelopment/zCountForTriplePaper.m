% This presents detailed methods that we need to have a record of but which does not belong in the paper:

% Method for counting nucleotides, basepairs, triples, etc.

NRList = 'Nonredundant_LT_4A_2011-01-07_list';
Verbose = 1;

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(NRList,0,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(NRList,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

Verbose = 0;
UsingLibrary = 1;
GUIActive = 0;
UsingLibrary = 0;

clear SearchFile

for i = 1:length(File),
  SearchFile{i} = File(i).Filename;
end

% -------------------------------------------------- Size of list

fprintf('There are %d 3D structures in %s\n', length(File), NRList);

%   1. Number of nucleotides.  FR3D search of NR 4A dataset from 2010-05-19.  Symbolic search, 2 nucleotides, no base constraint, distance constraint "=1 >" so that it finds all pairs of nucleotides one unit apart, exclude overlaps which excludes instances that are in redundant chains in the same PDB file.  Each pair found adds one distinct nucleotide.  This search finds 29907 pairs.  However, the last nucleotide of each file is not found, so we need to add the number of files in the NR dataset, which is 293.  Adding these gives 30,200.

% This method is unreliable, since it sometimes does and sometimes does not count the last nucleotide in a chain, depending on whether it is within 30 Angstroms of the first nucleotide in the next chain.

if 0 > 1,
  clear Query Search
  Query.ExcludeOverlap = 1;
  Query.SearchFiles = SearchFile;
  Query.Diff{1,2} = '=1 <';
  Query.Name = 'Nucleotides in NR set';
  Query = xConstructQuery(Query);
  Filenames = SearchFile;
  xFR3DSearch

  [L,N] = size(Search.Candidates);

  fprintf('Found %6d distinct nucleotides in %s\n', L+length(File), NRList);
end

% This is a more reliable method to count nucleotides in the NR dataset.  It puts all nucleotides from each file into a Search, then runs the program to exclude nucleotides from redundant chains.

clear Query Search
clear Candidates Discrepancy
Candidates = [];

for f = 1:length(File),
  NumNT = length(File(f).NT);
  Candidates = [Candidates; [(1:NumNT)' f*ones(NumNT,1)]];
end

Discrepancy = ones(size(Candidates(:,1)));
[Candidates] = xExcludeRedundantCandidates(File,Candidates,Discrepancy);

NumNRNT = length(Candidates(:,1));
fprintf('Found %6d distinct nucleotides in %s\n', length(Candidates(:,1)), NRList);

%---------------------------------------------------------------------------

%   2. Number of nucleotides involved in at least one basepair.  FR3D search of NR 4A dataset from 2010-05-19.  Symbolic search, 2 nucleotides, no distance constraint, basepair constraint "pair" (which gives basepair categories 1 through 12, and which excludes bifurcated basepairs).  This gives 25590 instances.  However, some nucleotides make more than one basepair, so this is not the final answer.  Next, I loaded the SearchSaveFile with Matlab and found the number of distinct nucleotide 1 and file number combinations.  This number is 

if 0 > 1,

% This method gives more candidates, not quite sure why, but it must have
% to do with redundant chains.

clear Query Search
Query.ExcludeOverlap = 1;
Query.SearchFiles = SearchFile;
Query.Edges{1,2} = 'pair';
Query.Diff{1,2} = '<';
Query.Name = 'Ordered basepairs in NR set';
Query = xConstructQuery(Query);
Filenames = SearchFile;
xFR3DSearch

[L,N] = size(Search.Candidates);

fprintf('Found %6d basepairs in %s\n', L, NRList);

AllPairNT = [Candidates(:,[1 3]); Candidates(:,[2 3])];
DistinctNT = unique(AllPairNT,'rows');

fprintf('Found %6d nucleotides (%7.4f%%) making at least one basepair in %s\n', length(DistinctNT(:,1)), 100*length(DistinctNT(:,1))/NumNRNT, NRList);

S2 = xDistinctCandidates(Search,2);
end


clear Query Search
Query.ExcludeOverlap = 1;
Query.SearchFiles = SearchFile;
Query.Edges{1,2} = 'pair';
Query.Name = 'Basepairs in NR set';
Query = xConstructQuery(Query);
Filenames = SearchFile;
xFR3DSearch

[S2,NumPairs] = xDistinctCandidates(Search);

fprintf('Found %6d basepairs in %s\n', NumPairs, NRList);

AllPairNT = [Candidates(:,[1 3]); Candidates(:,[2 3])];
DistinctNT = unique(AllPairNT,'rows');

fprintf('Found %6d nucleotides (%7.4f%%) making at least one basepair in %s\n', length(DistinctNT(:,1)), 100*length(DistinctNT(:,1))/NumNRNT, NRList);


%---------------------------------------------------------------------------
%   3. Number of distinct base triples with at least two basepairs.  FR3D search of NR 4A dataset from 2010-05-19.  Symbolic search, 3 nucleotides

clear Query Search
Query.ExcludeOverlap = 1;
Query.SearchFiles = SearchFile;
Query.Edges{1,2} = 'pair';
Query.Edges{2,3} = 'pair';
Query.Edges{1,3} = '~stack';
Query.Name = 'Base triples in NR set';
Query = xConstructQuery(Query);
Filenames = SearchFile;
xFR3DSearch

Search = xDistinctCandidates(Search);       % 

fprintf('Found %6d distinct base triples with at least two basepairs in %s\n', length(Search.Candidates(:,1)), NRList);

AllTripleNT = [Candidates(:,[1 4]); Candidates(:,[2 4]); Candidates(:,[3 4])];
DistinctNT = unique(AllTripleNT,'rows');

fprintf('Found %6d nucleotides (%7.4f%%) involved in at least one base triple in %s\n', length(DistinctNT(:,1)), 100*length(DistinctNT(:,1))/NumNRNT, NRList);

DistinctNT = unique(Candidates(:,[2 4]),'rows');

fprintf('Found %6d nucleotides (%7.4f%%) which make two or more basepairs simultaneously in %s\n', length(DistinctNT(:,1)), 100*length(DistinctNT(:,1))/NumNRNT, NRList);

% ------------------------------------ Look up the classifications

[L,N] = size(Search.Candidates);
N = N - 1;

for c = 1:L,                                % record interactions in triples
  f = Search.Candidates(c,N+1);             % file number
  Indices = Search.Candidates(c,1:N);
  E = File(f).Edge(Indices,Indices);
  Search.Edge{c} = E;
  Search.Data(c,:) = [File(f).NT(Indices).Code E(1,2) E(1,3) E(2,3)];
end

numbp = [];

for c = 1:L,
  e = abs(Search.Data(c,4:6));               % interactions
  numbp(c) = sum((e < 13) .* (e > 0));
  if numbp(c) < 2,
    fprintf('Wrong number of basepairs in candidate %d\n', c);
  end
end

fprintf('Found %6d distinct base triples with exactly two basepairs in %s\n', length(find(numbp==2)), NRList);

% Number of distinct base triples with three basepairs.

fprintf('Found %6d distinct base triples with three basepairs in %s\n', length(find(numbp==3)), NRList);

% ------------------------------------------------ Making three or more pairs

clear Query Search
Query.ExcludeOverlap = 1;
Query.SearchFiles = SearchFile;
Query.Edges{1,2} = 'pair';
Query.Edges{1,3} = 'pair';
Query.Edges{1,4} = 'pair';
Query.Edges{2,3} = '~stack';
Query.Edges{2,4} = '~stack';
Query.Edges{3,4} = '~stack';
Query.Name = 'Base making three pairs';
Query = xConstructQuery(Query);
Filenames = SearchFile;
xFR3DSearch

DistinctNT = unique(Candidates(:,[1 5]),'rows');

fprintf('Found %6d nucleotides (%7.4f%%) which make three or more basepairs simultaneously in %s\n', length(DistinctNT(:,1)), 100*length(DistinctNT(:,1))/NumNRNT, NRList);




% ------------------------------------------------ Coplanar search

clear Query Search
Query.ExcludeOverlap = 1;
Query.SearchFiles = SearchFile;
Query.Edges{1,2} = 'cp ncp';
Query.Edges{2,3} = 'cp ncp';
Query.Edges{1,3} = '~stack';
Query.Name = 'Coplanar';
Query = xConstructQuery(Query);

[File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
File = File(SIndex);

xFR3DSearch
Search = xDistinctCandidates(Search);       % 



% ------------------------------------ Look up the classifications

[L,N] = size(Search.Candidates);
N = N - 1;

fprintf('Found %d co-planar triples in %s\n', L, NRList);

for c = 1:L,                                % record interactions in triples
  f = Search.Candidates(c,N+1);             % file number
  Indices = Search.Candidates(c,1:N);
  E = File(f).Edge(Indices,Indices);
  Search.Edge{c} = E;
  Search.Data(c,:) = [File(f).NT(Indices).Code E(1,2) E(1,3) E(2,3)];
end

% full(sortrows(Search.Data,[1 2 3 4 5 6 1 2 3]))

k = find((Search.Data(:,4)-fix(Search.Data(:,4))) == 0.5);

%fprintf('Subcategory e basepairs in triples\n');
%full(sortrows(Search.Data(k,:),[4 5 6 1 2 3]))

numbp = [];

for c = 1:L,
  e = abs(Search.Data(c,4:6));               % interactions
  numbp(c) = sum((e < 13) .* (e > 0));
end

fprintf('Found %6d distinct co-planar triples with exactly two basepairs in %s\n', length(find(numbp==2)), NRList);

% Number of distinct base triples with three basepairs.

fprintf('Found %6d distinct co-planar triples with three basepairs in %s\n', length(find(numbp==3)), NRList);






break






% --------------------------------------------------------------------------
% Coplanar Triples from original search
% Evaluate updated basepair classification

load 2010-07-20_11_45_16-cp_ncp_not_stack_triples_excl_1QCU.mat

Filenames = Search.Filenames;
[File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
File = File(SIndex);

% ------------------------------------ Look up the classifications

[L,N] = size(Search.Candidates);
N = N - 1;

for c = 1:L,                                % record interactions in triples
  f = Search.Candidates(c,N+1);             % file number
  Indices = Search.Candidates(c,1:N);
  if max(Indices) <= File(f).NumNT,
    E = File(f).Edge(Indices,Indices);
    Search.Edge{c} = E;
    Search.Data(c,:) = [File(f).NT(Indices).Code E(1,2) E(1,3) E(2,3)];
  else
    fprintf('Size of %s has changed\n',File(f).Filename);
    E = Search.File(f).Edge(Indices,Indices);
    Search.Edge{c} = E;
    Search.Data(c,:) = [Search.File(f).NT(Indices).Code E(1,2) E(1,3) E(2,3)];
  end
end

k = find( mod(abs(Search.Data(:,4)),100) < mod(abs(Search.Data(:,6)),100));
length(k)

numbp = [];

for c = 1:L,
  e = abs(Search.Data(c,4:6));               % interactions
  numbp(c) = sum((e < 13) .* (e > 0));
end

fprintf('Found %6d distinct co-planar triples with exactly two basepairs in %s\n', length(find(numbp==2)), NRList);

% Number of distinct base triples with three basepairs.

fprintf('Found %6d distinct co-planar triples with three basepairs in %s\n', length(find(numbp==3)), NRList);


break

% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% Triples in bacterial ribosomes

Verbose = 1;
UsingLibrary = 1;
GUIActive = 0;
UsingLibrary = 0;

clear Query Search
Filenames = {'3i8i','2qbe','1j5e','2qan'};
Query.ExcludeOverlap = 1;
Query.SearchFiles = Filenames;
Query.Edges{1,2} = 'pair';
Query.Edges{2,3} = 'pair';
Query.Edges{1,3} = '~stack';
Query.Name = 'Base triples in bacterial ribosomes;'
Query = xConstructQuery(Query);

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);
xFR3DSearch

Search = xDistinctCandidates(Search);       % 

% ------------------------------------ Look up the classifications

[L,N] = size(Search.Candidates);
N = N - 1;

for c = 1:L,                                % record interactions in triples
  f = Search.Candidates(c,N+1);             % file number
  Indices = Search.Candidates(c,1:N);
  E = File(f).Edge(Indices,Indices);
  Search.Edge{c} = E;
  Search.Data(c,:) = [File(f).NT(Indices).Code E(1,2) E(1,3) E(2,3)];
  Search.Chain(c,:) = cat(2,File(f).NT(Indices).Chain);
end

numbp = [];

for c = 1:L,
  e = abs(Search.Data(c,4:6));               % interactions
  numbp(c) = sum((e < 13) .* (e > 0));
  if numbp(c) < 2,
    fprintf('Wrong number of basepairs in candidate %d\n', c);
  end
end

for f = 1:length(Filenames),
  C = cat(2,File(f).NT.Chain);
  U = unique(C);

  for u = 1:length(U),
    k = find((Search.Candidates(:,N+1) == f) .* (Search.Chain(:,1) == U(u)));
    kk = find((Search.Candidates(:,N+1) == f) .* (numbp' == 2) .* (Search.Chain(:,1) == U(u)));
   
    fprintf('Found %4d base triples in %s chain %s, of which %4d have 2 basepairs and %4d have 3.\n', length(k), Filenames{f}, U(u), length(kk), length(k)-length(kk));

  end
end

% ------------------------------------------- Write to spreadsheet

clear T

for c = 1:L,
  f = Search.Candidates(c,N+1);
  T{c,1} = Filenames{f};
  ci = sort(Search.Candidates(c,1:3));
  T{c,12} = '';
  for n = 1:3,
    T{c,2*n}   = File(f).NT(ci(n)).Base;
    T{c,2*n+1} = File(f).NT(ci(n)).Number;
    T{c,12}    = [T{c,12} File(f).NT(ci(n)).Chain];
  end
  T{c,8}   = zEdgeText(File(f).Edge(ci(1),ci(2)));
  T{c,9}   = zEdgeText(File(f).Edge(ci(1),ci(3)));
  T{c,10}  = zEdgeText(File(f).Edge(ci(2),ci(3)));
  if length(File(f).NT) > 2000,
    T{c,11} = '5S_23S';
  else
    T{c,11} = '16S';
  end
end 

T

break

%T = sortrows(T,[11 3 5 7]);
xlswrite('Triples in ribosomes.xls',T)
%xlswrite('Triples with near pairs in ribosomes.xls',T)

