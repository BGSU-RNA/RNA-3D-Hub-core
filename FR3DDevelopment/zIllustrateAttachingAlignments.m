
% This is an illustration of how to attach an alignment to a 3D structure

File = zAddNTData('2AVY');

File = zAttachAlignment(File,1,{'16S-simplified.fasta'},{'A'},354);

% here is another, with Anton's motifs

load MotifLibrary\IL_032.mat

for f = 1:length(Search.File),
  Filenames{f} = Search.File(f).Filename(1:4);
end

FN = unique(Filenames);

File = zAddNTData(FN);

[L,N] = size(Search.Candidates);
N = N - 1;

for c = 1:L,
  ff = Search.Candidates(c,N+1);               % original file number
  f  = find(ismember(FN,Filenames{ff}));       % new file number
  NewCand(c,N+1) = f;
  for n = 1:N,                                 % loop through nucleotides
    i = Search.Candidates(c,n);
    NewCand(c,n) = zIndexLookup(File(f),Search.File(ff).NT(i).Number,Search.File(ff).NT(i).Chain);
  end
end

S.Candidates = NewCand;
S.Filenames = FN;

xDisplayCandidates(File,S);

Verbose = 1;

disp('Make the working directory be Alignments');

File = zAttachAlignment(File,1);
[File,seqs,codes,counts] = zGetVariantsFromAlignments(File,S,Verbose)
