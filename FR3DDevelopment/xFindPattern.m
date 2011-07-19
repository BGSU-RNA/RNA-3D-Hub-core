% xFindPattern hijacks FR3D to search for a pattern in a field of points

% function [Candidates] = xFindPattern(field,pattern,D,Verbose)

if nargin < 4,
  Verbose = 2;
end

Verbose = 0;

parameters = 1:7;

for p = 1:length(parameters),

parameter = parameters(p);

switch parameter,

case 1,
  N = 4;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 20;                               % number of extra points in the field
  D = 0.02;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 2,
  N = 4;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 200;                              % number of extra points in the field
  D = 0.05;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 3,
  N = 4;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 2000;                             % number of extra points in the field
  D = 0.05;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 4,
  N = 5;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 20;                               % number of extra points in the field
  D = 0.02;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 5,
  N = 5;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 200;                              % number of extra points in the field
  D = 0.05;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 6,
  N = 5;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 2000;                             % number of extra points in the field
  D = 0.02;                             % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

case 7,
  N = 5;                                % number of points in pattern
  d = 3;                                % dimensionality of the data
  F = 2000;                             % number of extra points in the field
  D = 0.1;                              % maximum discrepancy
  s = RandStream.create('mt19937ar','seed',5489);

end

settings(parameter,:) = [N d F D];

RandStream.setDefaultStream(s);

pattern = rand(N,d);                  % pattern to look for
field = rand(F,d)*(F^(1/d));          % where to look; keep density constant
field = [field; pattern];             % make sure the pattern is in the field
pattern = round(pattern*1000)/1000;   % round to three digits
field = round(field*1000)/1000;       % round to three digits

if Verbose > 0,
  field
  pattern
end

colormap('default');
map = colormap;

if Verbose > 1,
  figure(1)
  clf
  subplot(3,1,1);
  plot(field(:,1),field(:,2),'.k');
  axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
  ylabel('Field');
  subplot(3,1,2);
  scatter(pattern(:,1),pattern(:,2),30,1:N,'filled');
  colormap(map(8:56,:));
  axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
  ylabel('Pattern');
end

% -------------------------------- set up a fictitious file to search in

clear File
File.Distance = zDistance(field);
File.Filename = 'xFindPattern';
File.NumNT    = length(field(:,1));
for n = 1:File.NumNT,
  File.NT(n).Code = 1;
end

% -------------------------------- set up the pattern as the query for FR3D

clear Query
Query.Geometric = 1;
Query.Diameter  = 8;
Query.Distance  = zDistance(pattern);
Query.NumNT     = length(pattern(:,1));
Query.LocWeight = ones(1,Query.NumNT);
Query.AngleWeight = ones(1,Query.NumNT);
Query.DiscCutoff = D;
Query.RelCutoff = Query.DiscCutoff;
Query.SSCutoff  =(Query.NumNT^2)*(Query.RelCutoff^2)*cumsum(Query.LocWeight);

% ------- The next line makes sure that no subset screening will occur:

Query.SSCutoff(4:end) = Inf*Query.SSCutoff(4:end);

% -------------------------------------- Find candidates
t = cputime;

[Candidates,Lists] = xFindCandidates(File,Query,12);

timetest(parameter) = cputime - t;     % store the elapsed time

fprintf('Pattern %d has %d points in field of %d discrep %0.2f\n',parameter,N,F+N,D);

fprintf('Found %d candidates.\n', length(Candidates(:,1)));

% -------------------------------------- Assemble pairwise screens into lists

ListMatrix = [];
for p = 1:N,
  for q = (p+1):N,
    [i,j] = find(Lists{p,q});
    NLM = [p*ones(size(i)) q*ones(size(i)) i j];
    ListMatrix = [ListMatrix; NLM];
  end
end

Candidates = sortrows(Candidates);
ListMatrix = sortrows(ListMatrix);

% ----------------------------------- Write out to data files

FNBase = sprintf('%d_points_in_field_of_%d_discrep_%0.2f',N,F+N,D);
FNBase = strrep(FNBase,'.',',');

zWriteMatrixAsCSV(pattern,[FNBase '_pattern.txt'],1);
zWriteMatrixAsCSV(field,[FNBase '_field.txt'],1);
zWriteMatrixAsCSV(Candidates(:,1:N),[FNBase '_candidates.txt'],0);
zWriteMatrixAsCSV(ListMatrix,[FNBase '_lists.txt'],0);

% ----------------------------------------- Sort by discrepancy

if Verbose > 1,

  L = length(Candidates(:,1));                     % number of candidates

  fprintf('Calculating discrepancies: ');

  MCC = pattern - ones(Query.NumNT,1)*mean(pattern);

  Discrepancy = [];

  for i=1:L,
    C  = field(Candidates(i,1:Query.NumNT),:);         % coordinates of candidate i
    CC = C - ones(Query.NumNT,1)*mean(C);

    R = zBestRotation(CC, diag(Query.LocWeight)*MCC);      % candidate onto model
  
    S = Query.LocWeight * sum(((MCC - CC*R').^2)')';  % distances between centers

    Discrepancy(i) = S;

    if (mod(i,round(L/10)) == 0) && (Verbose > 0)
      fprintf(' %d', fix((L-i)*toc/i)); 
      drawnow
    end
  end

  fprintf('\n');

  Discrepancy = sqrt(Discrepancy)/Query.NumNT;

  [y,i]       = sort(Discrepancy);                    % sort by discrepancy
  Candidates  = Candidates(i,:);
  Discrepancy = Discrepancy(i);

  fprintf('Calculating discrepancy took        %8.3f seconds\n',toc);

  for j = 1:L,
    figure(1)
    subplot(3,1,3);
    i = Candidates(j,1:N);
    scatter(field(i,1),field(i,2),30,1:N,'filled');
    axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
    colormap(map(8:56,:));
    ylabel('Candidate');
    fprintf('Candidate %d of %d has discrepancy %7.4f\n', j, L, Discrepancy(j));
    pause
  end
end

end

if length(parameters) > 1,
  ttid = fopen('xFindPattern_runtimes.txt','w');
  for p = 1:length(parameters),
    parameter = parameters(p);
    fprintf(ttid,'Pattern %d has %d points in field of %d discrep %0.2f.  Runtime %12.8f seconds\n',parameter,settings(parameter,1),settings(parameter,1)+settings(parameter,3),settings(parameter,4),timetest(parameter));
  end
  fclose(ttid);
end
