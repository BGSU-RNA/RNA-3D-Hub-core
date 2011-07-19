% xListPairsByIDI takes the results of a search for 2-NT motifs and sorts them by IDI

Verbose = 1;

[IDI,Candidates] = xRankCandidatesByIDI(Search.File,Search.Candidates,Verbose);

Search.Candidates = Candidates;
Search.Discrepancy = IDI;

xListCandidates(Search)
xDisplayCandidates(Search.File,Search)
