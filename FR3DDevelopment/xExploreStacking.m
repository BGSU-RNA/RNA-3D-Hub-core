% given a motif, do searches on all the stackings and see what sequence variability they exhibit when ranked by IDI

Verbose = 0;

load LIB00003_IL_tSH-tWH-tHS.mat

load 2008-12-12_17_28_25-stacked_cww.mat

S = Search;

[L,N] = size(S.Candidates);
N = N - 1;

f = S.Candidates(1,N+1);                   % File of search
F = S.File(f);

%Filenames = 'HighResolution_list';
Filenames = 'NonRedundant_2008_02_21_list';
% File = zAddNTData(Filenames);

fprintf('The motif being studied is from %s.\n', S.File(f).Filename);

IDIcutoff = 2;

clear UsingLibrary
GUIactive = 1;

Letters = 'ACGU';

for ii = 1:(N-1),
  for jj = (ii+1):N,
    e = abs(F.Edge(S.Candidates(1,ii),S.Candidates(1,jj)));

    if e > 0 && e < 24,
      Query.Description    = 'Pairwise interaction in motif';
      Query.Name           = 'Pairwise interaction in motif';
      Query.Filename       = F.Filename;
      Query.NTList         = {F.NT(S.Candidates(1,ii)).Number F.NT(S.Candidates(1,jj)).Number};
      Query.ChainList      = {F.NT(S.Candidates(1,ii)).Chain F.NT(S.Candidates(1,jj)).Chain};
      Query.DiscCutoff     = 0.7;
      Query.SearchFiles    = {Filenames};

      xFR3DSearch

      [IDI,Candidates] = xRankCandidatesByIDI(Search.File,Search.Candidates,Verbose);

      Search.Candidates = Candidates;
      Search.Discrepancy = IDI;

      counts = zeros(4,4);

      c = 1;
      while IDI(c) < IDIcutoff,
        f  = Candidates(c,3);
        i1 = Candidates(c,1);
        i2 = Candidates(c,2);
        c1 = Search.File(f).NT(i1).Code;
        c2 = Search.File(f).NT(i2).Code;
        counts(c1,c2) = counts(c1,c2) + 1;
        c = c + 1;
      end

      f  = Candidates(1,3);
      i1 = Candidates(1,1);
      i2 = Candidates(1,2);
      N1 = Search.File(f).NT(i1);
      N2 = Search.File(f).NT(i2);

      if e < 14,
        te = 'basepaired';
      else
        te = 'stacked';
      end

      fprintf('Among %s nucleotides within IDI %4.2f of %s%s - %s%s %s from the query motif, here are the counts of base combinations:\n', te, IDIcutoff, N1.Base,N1.Number,N2.Base,N2.Number,zEdgeText(e,0,N1.Code,N2.Code));

      fprintf('        A      C      G      U  (corresponding to %s%s)\n',N2.Base,N2.Number);
      for i = 1:4,
        fprintf('%c   %5d  %5d  %5d  %5d\n', Letters(i), counts(i,1), counts(i,2), counts(i,3), counts(i,4));
      end
      fprintf('\n');

%      xListCandidates(Search,100);
%      xDisplayCandidates(Search.File,Search);
    end
  end
end


