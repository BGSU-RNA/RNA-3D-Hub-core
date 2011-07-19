% This shows how to list candidates from a search by similarity, with the candidates annotated according to what motif they're part of

% in Matlab, do the following:

% load 2009-03-12_17_53_43-Internal_loops_avoiding_ncWW.mat
% load 2009-03-12_20_44_20-Hairpin_loops_with_extra_bases.mat

File = zAddNTData(Search.Filenames);

File = xAnnotateWithKnownMotifs(File,1,0);

Limit = length(Search.Candidates(:,1));
L     = length(Search.Candidates(:,1));

Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies

clear p

p(1:Limit) = zOrderbySimilarity(Search.Disc(1:Limit,1:Limit));

Search2             = Search;
Search2.Candidates  = Search.Candidates(p,:);
Search2.Discrepancy = Search.Discrepancy(p);
Search2.Marked      = Search.Marked(p);
Search2.Disc        = Search.Disc(p,p);
Search2.DiscComputed= Search.DiscComputed(1,p);
Search2.File        = File;                       % use the full files

xListCandidates(Search2,Inf);
