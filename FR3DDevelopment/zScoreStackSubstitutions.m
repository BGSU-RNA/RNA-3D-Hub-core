% zScoreStackSubstitutions(File,Index1,Index2,Verbose,Stacks) calculates a 4x4 substitution score matrix for the stack between nucleotides Index1 and Index2 in File.  You need to make sure that Index1 and Index2 are stacked.  




% Run it this way:

% [M,Stack] = zScoreStackSubstitutions('1S72','1106','1107',2);
% [M,Stack] = zScoreStackSubstitutions('1S72','1106','1107',2,Stack);


% Stacks = [];
% for i = 1:allstacks,
%   [Subs,Stacks] = zScoreStackSubstitutions(File,Index1,Index2,Verbose,Stacks);
% end

% The first time, it will load the Stacks variable and pass it back.
% This saves time.


function [Subs,Stacks] = zScoreStackSubstitutions(File,Index1,Index2,Verbose,Stacks)

if nargin < 5,
  load 2010-05-30_10_07_51-s35_NR_4A.mat
  Stacks{1} = Search;
  load 2010-05-29_10_24_40-s33_NR_4A_exclude_redundant.mat
  Stacks{2} = Search;
  load 2010-05-29_10_39_07-s55_NR_4A.mat
  Stacks{3} = Search;
end

if nargin < 4,
  Verbose = 0;
end

if isempty(Stacks),
  load 2010-05-30_10_07_51-s35_NR_4A.mat
  Stacks{1} = Search;
  load 2010-05-29_10_24_40-s33_NR_4A_exclude_redundant.mat
  Stacks{2} = Search;
  load 2010-05-29_10_39_07-s55_NR_4A.mat
  Stacks{3} = Search;
end

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if Index1 is a string, look up the index

if strcmp(class(Index1),'char'),
  Index1 = zIndexLookup(File,Index1);
end

if strcmp(class(Index2),'char'),
  Index2 = zIndexLookup(File,Index2);
end

% ---------------------------------------------- Identify type of stack

Trouble = 0;

      switch fix(File.Edge(Index1,Index2)),
      case 21,
        S = Stacks{1};                          % s35
      case {22, -22},
        S = Stacks{2};                          % s33
      case {23, -23},
        S = Stacks{3};                          % s55
      case -21,
        S = Stacks{1};                          % s53
        S.Candidates = S.Candidates(:,[2 1 3]);   % swap nucleotide order
      case 121,																	% Added by Steve to deal with near stacking interactions
        S = Stacks{1};                          % s35							% Added by Steve to deal with near stacking interactions
      case {122, -122},															% Added by Steve to deal with near stacking interactions
        S = Stacks{2};                          % s33							% Added by Steve to deal with near stacking interactions
      case {123, -123},															% Added by Steve to deal with near stacking interactions
        S = Stacks{3};                          % s55							% Added by Steve to deal with near stacking interactions
      case -121,																% Added by Steve to deal with near stacking interactions
        S = Stacks{1};                          % s53							% Added by Steve to deal with near stacking interactions
        S.Candidates = S.Candidates(:,[2 1 3]);   % swap nucleotide order		% Added by Steve to deal with near stacking interactions
      otherwise,
        disp('What case is this?');
        fix(File.Edge(Index1,Index2))  % Steve changed this from Edge() to File.Edge()
        Trouble = 1;
      end

% ----------------------------------------------- Analyze database of stacks

if Trouble == 0,
      NS = length(S.Candidates(:,1));           % number of stacks
      for c = 1:NS,
        f = S.Candidates(c,3);                  % file number
        SC(c,1) = S.File(f).NT(S.Candidates(c,1)).Code;  % code of first base
        SC(c,2) = S.File(f).NT(S.Candidates(c,2)).Code;  % code of second base
      end

      H = zeros(4,4);

        clear Model

        Model.NT(1) = File.NT(Index1);
        Model.NT(2) = File.NT(Index2);

% zDisplayNT(Model)

% ----------------------------------------------- Compare to database of stacks

        [Discrepancy, Candidates, i] = xRankCandidatesIDI(S.File,Model,S.Candidates,0);

        NewSC = SC(i,:);                        % re-order candidate codes

        MinIDI = 40*ones(4,4);
        
        for c = NS:-1:1,
          MinIDI(NewSC(c,1),NewSC(c,2)) = Discrepancy(c);
        end

        G = 1 ./ (1+MinIDI.^2);
        G = G / sum(sum(G));

        Subs = G;

        if Verbose > 0,
          fprintf('We get this min IDI matrix for %s %s%s %s%s %s %s\n',File.Filename,Model.NT(1).Base,Model.NT(1).Number,Model.NT(2).Base,Model.NT(2).Number,zEdgeText(fix(File.Edge(Index1,Index2))));
 
MinIDI

          fprintf('We turn it into this matrix of scores:\n');
          G
        end

else

  Subs = NaN * ones(4,4);

end
