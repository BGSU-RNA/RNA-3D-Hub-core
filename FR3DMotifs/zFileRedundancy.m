% zFileRedundancy(Filenames) explores possible redundancy between PDB files listed in Filenames.  It sorts the files by number of nucleotides, then compares files with similar numbers of nucleotides, performing a Needleman-Wunsch alignment of their bases.  If the percentage of bases which align exceeds the parameter p, the pair is kept for further examination.  It groups together structures connected by chains of greater than p percent sequence identify, then geometrically superimposes all pairs in each group, then prints a report so a human can decide which structures to keep.

% zFileRedundancy('Allfiles_list') % should give a huge report
% zFileRedundancy('NonRedundant_2008_02_21_list') % should show very little possible redundancy

function [t,n] = zFileRedundancy(reportdate,t,n)

if nargin < 3,
     load PDBInfo
%    load([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat']);
end

% Columns of t correspond to this list:
% 1 A	Structure ID
% 2 B	Descriptor / Structure title
% 3 C	Experimental Technique / Exp. method
% 4 D	Release Date
% 5 E	Authors
% 6 F	Keywords
% 7 G	Resolution (Å)
% 8 H	Source

% We should add an additional column for NDB ID
% 11    Sequence of longest chain
% 12    Best chains, as determined by zBestChains

% In zFileRedundancy, column 10 is set
% 10    PDB ID of structure which represents this one

% Columns of n correspond to this list:
% n(i,1)    Resolution (NaN if NMR)
% n(i,2)    Number of nucleotides
% n(i,3)    Number of basepairs including bifurcated
% n(i,4)    Number of nucleotides in longest chain
% n(i,5)    Starting index of longest chain
% n(i,6)    Ending index of longest chain
% n(i,7)    Number of nucleotides in NR chains
% n(i,8)    Number of basepairs in NR chains
% n(i,9)    Number of cWW basepairs in structure
% n(i,10)   Number of non-cWW basepairs in structure


if ~exist('reportdate')
  reportdate = datestr(date,'yyyy-mm-dd');
end

tot = cputime;                    % keep track of total time

p = 0.95;                         % cutoff base match fraction

maxd = 0.5;                       % cutoff discrepancy between structures
NTLimit = 30000;                  % above this limit, do not align sequences
MaxRes  = 100;                    % maximum resolution value to use
SL = 4000;                        % upper limit on # bases to align

Criterion = 5;
                                  % 1-earliest release date
                                  % 2-resolution
                                  % 3-number of nucleotides
                                  % 4-#pairs
                                  % 5-highest #BP / #nucleotides
                                  % 6-most recent release date

                                  % add 10, use preferred list to override

OriginalNumStruct = length(t(:,1));

diary(['Redundancy_Report_' reportdate '_Criterion_' num2str(Criterion) '.txt']);

Preferred = {'1J5E'};
Preferred = {};

JoinList = {'1FG0','1S72'};       % manually declare these to be redundant
%JoinList = [JoinList; {'1C2W','2AW4'}]; % more can be added

% --------------------------------- Look up information on these structures

fprintf('Starting with %4d structures from the PDB.\n', length(t(:,1)));

[y,i] = sort(n(:,4));         % sort by number of nucleotides in longest chain
t = t(i,:);
n = n(i,:);

F = length(t(:,1));                    % number of files

fprintf('Preparing a redundancy report on %d RNA 3D structures.\n',F);
fprintf('Sequence identity cutoff will be %7.2f.\n', p);
fprintf('Geometric discrepancy cutoff will be %7.2f.\n', maxd);

% --------------------------------- Compare sequences between files

Close = 0.5*sparse(eye(F));       % indicator of whether sequences are close
prop  = 0.5*sparse(eye(F));

clear align
align{F,F} = [];

stim = cputime;
tim  = cputime;

Linked = zeros(1,F);              % whether each file is already redundant

for i = 1:(F-1),                          % loop through all files
 if n(i,4) < 2,                           % such files are too short to align
                                          % do nothing
 elseif Linked(i) <= 2 || n(i,4) < 1300,      % if already linked 0, 1, 2 times
                                          % or if a small structure

  if n(i,4) < 19,
    maxlength = n(i,4);                   % short chains need to be identical in length
  else
    maxlength = 2*n(i,4);                 % otherwise up to twice the length
  end

  J = find( ((1:F)' > i) .* (n(1:F,4) <= maxlength) .* n(1:F,4) >= n(i,4));
                                          % others w/chain length in right range

  if (cputime - tim > 10) || (mod(i,20) == 0) || (i==1),
    fprintf('Comparing sequence from structure %4d, %s which has %4d nucleotides in its longest chain, to %3d others of a similar size\n', i, t{i,1}, n(i,4), length(J));
    tim = cputime;
  end

  for k = 1:length(J),                        % loop through other structures
    j = J(k);                                 % look at structure j
    if n(i,2) <= NTLimit || length(t{i,8}) == 0 || length(t{j,8}) == 0,

      in = lower(t{i,8});                     % biological source
      jn = lower(t{j,8});                     % biological source

      CompareOrg = zCompareOrganismNames(in,jn);  % 1 if definitely diff't,
                                              % 0 if could be the same

      ni = n(i,10);                           % number of non-cWW pairs
      nj = n(j,10);                           % number of non-cWW pairs

      if (CompareOrg > 0) || ((min(ni,nj) == 0) && max(ni,nj) > 0),
        pro = 0;                              % do not link these files
      else
        ti = t{i,11};                           % seq from longest chain
        ti = ti(1:min(length(ti),SL));          % only compare first SL bases
        tj = t{j,11};                           % seq from longest chain
        tj = tj(1:min(length(tj),SL));          % only compare first SL bases
        [matches,a,b,ss,tt] = zNeedlemanWunsch(ti,tj); % sequence alignment
        e = (t{i,11}(a) == t{j,11}(b));         % locations of agreement
        pro = matches/min(length(ti),length(tj));    % percent identity

        % -------------------------------- Record close sequence agreement

        if (length(ti) < 10 && matches == length(ti)) || ...
           (length(ti) >= 10 && min(length(ti),length(tj)) - matches <= 4) || ...
           (length(ti) >= 80 && pro >= p),
          Close(i,j) = 1;                       % sequences agree well enough
          prop(i,j)  = pro;
          Linked(j)  = Linked(j) + 1;           % tally how many times linked
          if n(i,2) > NTLimit,
            fprintf('Sequence identity with %s is %6.4f\n', t{j,1}, pro);
          end
          align{i,j} = [a(e); b(e)];        % save alignment data for later
        end

        if prop(i,j) > 0.9,                 % if even a little sequence agreement
          align{i,j} = [a(e); b(e)];        % save alignment data for later
        end
      end

    end
    j = j + 1;
  end
 end

end

fprintf('\nComparing sequences from structures took %7.1f minutes.\n\n', (cputime-stim)/60);

% ------------------------------------------ Join specified structures

if exist('JoinList'),
  for j = 1:length(JoinList(:,1)),
    a = find(ismember(t(:,1),JoinList{j,1}));
    b = find(ismember(t(:,1),JoinList{j,2}));
    Close(a,b) = 1;
    prop(a,b) = 1;
  end
end

fprintf('\nExtending by symmetry and transitivity ');

ttime = cputime;

Close = Close + Close';           % extend by symmetry
Close = double(Close > 0);        % just make it 0-1 valued
prop  = prop + prop';

closeseq = Close;                 % store this for later

different = 1;
while different > 0,
  Previous = Close;
  Close = Close * Close;
  Close = double(Close > 0);            % make into 0 and 1
%  max(max(abs(Previous - Close)))
  if max(max(abs(Previous - Close))) == 0,
    different = 0;
  end
end

clear Previous

fprintf('took %7.1f minutes.\n\n', (cputime-ttime)/60);

figure(2)
clf
spy(Close)
title(['Structures connected by a chain of more than ' num2str(p*100) '% similarity']);
drawnow

save('FileRedundancy.mat')
whos

% ---------------------------------------- Superimpose geometrically

fprintf('Superimposing structures with sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered
Discrepancies = sparse(F,F);

stim = cputime;

for i = 1:(F-1),                  % loop through files
  if done(i) == 0,                % if file i hasn't been considered yet
    j = find(Close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;                % file i has a unique sequence
    else

      clear File

      for jj = 1:length(j),
        FF = zAddNTData(t{j(jj),1},0,[],1);
        FFF.NT = FF.NT;
        FFF.NumNT = FF.NumNT;
        FFF.Filename = FF.Filename;
        FFF.Info = FF.Info;
        FFF.LongestChain = FF.LongestChain;
        File(jj) = FFF;                         % only keep this data
      end
      clear FF FFF

      for m = 1:length(j),
        fprintf('%4s has %4d nucleotides, %4d in longest chain, %4d basepairs, ', File(m).Filename, File(m).NumNT, n(j(m),4), n(j(m),3));
        if isempty(File(m).Info.Resolution),
          fprintf('resolution  ---- ');
        else
          fprintf('resolution %5.2f ', File(m).Info.Resolution);
        end

        Info = File(m).Info;

        fprintf(' %10s | %s | %s | %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
      end

      fprintf('Percent agreement of base sequence, using alignment of sequences\n');

      fprintf('           ');
      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        done(j(m)) = 1;
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end
        for nn = 1:length(j),
          if nn == m || prop(j(m),j(nn)) == 0,
            T = [T sprintf('       ')];
          else
            T = [T sprintf(' %5.1f ', 100*full(prop(j(m),j(nn))))];
          end
        end
        fprintf('%s\n',T);
      end
      fprintf('\n');

      fprintf('Sequences of longest chains:\n');
      for m = 1:length(j),                        % print sequences
        if length(t{j(m),11}) < 10000,
          fprintf('%4s %s\n', File(m).Filename, t{j(m),11});
        end
      end
      fprintf('\n');


if length(j) > 30,
  fprintf('Calculating discrepancies\n');
  drawnow
end

      for m = 1:length(j),
        for nn = (m+1) : length(j),
          if ~isempty(align{j(m),j(nn)}),
            malign = align{j(m),j(nn)}(1,:);     % use stored data
            nalign = align{j(m),j(nn)}(2,:);
          else
            malign = [];
            nalign = [];
          end

          if isempty(malign) || closeseq(j(m),j(nn)) == 0,
            d = Inf;
          else
            try
              d = xDiscrepancy(File(m),File(m).LongestChain(1)-1+malign,File(nn),File(nn).LongestChain(1)-1+nalign);
            catch
              fprintf('Files %s and %s have trouble with their longest chains\n', t{j(m),1}, t{j(nn),1});
              d = Inf;
            end
          end

          Discrepancies(j(m),j(nn)) = d;
          Discrepancies(j(nn),j(m)) = d;

          if ~(d <= maxd),                  % allow for d = NaN
            closeseq(j(m),j(nn)) = 0;       % these are not that close!
            closeseq(j(nn),j(m)) = 0;
          end
        end
      end

      fprintf('Geometric discrepancy between aligned bases, where there is sequence similarity\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end


        for nn = 1:length(j),
          discr = full(Discrepancies(j(m),j(nn)));
          if nn == m,
            T = [T sprintf('       ')];
          elseif discr < Inf,
            T = [T sprintf(' %5.2f ', discr)];
          else
            T = [T sprintf('       ')];
          end
        end
        fprintf('%s\n',T);
      end
      fprintf('\n');
    end
  end
  drawnow
end

fprintf('\nSuperimposing structures took %7.1f minutes.\n\n', (cputime-stim)/60);

% Now repeat the analysis, having removed links between structures that are
% close in sequence but don't superimpose well geometrically.

fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n\n');

clear Text

Close = closeseq;                 % return to the previous relations

% ------------------------------------------ Join specified structures

if exist('JoinList'),
  for j = 1:length(JoinList(:,1)),
    a = find(ismember(t(:,1),JoinList{j,1}));
    b = find(ismember(t(:,1),JoinList{j,2}));
    Close(a,b) = 1;
    prop(a,b) = 1;
  end
end

different = 1;
while different > 0,
  Previous = Close;
  Close = Close * Close;
  Close = double(Close > 0);            % make into 0 and 1
  max(max(abs(Previous - Close)))
  if max(max(abs(Previous - Close))) == 0,
    different = 0;
  end
end

clear Previous

figure(3)
clf
spy(Close)
title(['Structures connected by ' num2str(p*100) '% sequence similarity and ' num2str(maxd) ' discrepancy']);
drawnow

for i = 1:F,
  fprintf('%4s is in a group with %8d files\n', t{i,1}, full(sum(Close(i,:) > 0)));
end

% -------------------------------------------

fprintf('Listing structures with more than %5.2f%% sequence similarity and less than %5.2f geometric discrepancy.\n',100*p,maxd);

fprintf('Structures are sorted by ');
switch mod(Criterion,10)
  case 1, fprintf('release date, earliest first\n');
  case 2, fprintf('resolution, best resolution first\n');
  case 3, fprintf('number of nucleotides (decreasing)\n');
  case 4, fprintf('number of basepairs (decreasing)\n');
  case 5, fprintf('number of basepairs per nucleotide (decreasing)\n');
  case 6, fprintf('release date, most recent first\n');
end

fprintf('The first structure listed will be chosen as the representative');
if Criterion > 10,
  fprintf(' unless one of these structures is found in the preferred list of structures.\n');
else
  fprintf('.\n');
end

fprintf('\n');

fprintf('Building equivalence classes for %d files\n', F);

done = zeros(1,F);                % whether each file has been considered
c = 0;                            % counts the number of files selected so far

for i = 1:F,                      % loop through all files
  if done(i) == 0,                % file does not already appear in a report
    j = find(Close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;                % no need to display this one again
      t{i,10} = t{i,1};           % this file represents itself
      j = 1;
    else

      crit = [];                   % criterion for sorting and choosing

      for f = 1:length(j),
        done(j(f)) = 1;            % no need to display this one again
                                   % sorting criteria
        crit(f,1) = datenum(t{j(f),4}, 'yyyy-mm-dd');
        crit(f,2) = n(j(f),1);            % resolution
        crit(f,3) = -n(j(f),2);           % number of nucleotides
        crit(f,4) = -n(j(f),3);               % number of pairs
        crit(f,5) = -n(j(f),3)/n(j(f),2); % ratio # bp / # nucleotides
%        crit(f,5) = round(crit(f,5)*50);  % round to nearest 0.02
      end

      fprintf('\n');

      switch mod(Criterion,10)
        case 1, [y,k] = sortrows(crit,[1 5 2 3 4]);
        case 2, [y,k] = sortrows(crit,[2 5]);
        case 3, [y,k] = sortrows(crit,[3 2]);
        case 4, [y,k] = sortrows(crit,[4 5]);
        case 5, [y,k] = sortrows(crit,[5 2 -1]);
      end

      j = j(k);                            % order by the criterion

      for m = 1:length(j),
        fprintf('%4s has %4d nucleotides, %4d pairs, ', t{j(m),1}, n(j(m),2), n(j(m),3));
        fprintf('resolution %5.2f ', n(j(m),1));
        fprintf(' %10s | %s | %s | %s\n', t{j(m),4}, t{j(m),8}, t{j(m),2}, t{j(m),5});
      end

      if n(i,2) <= NTLimit,

        fprintf('Percent agreement of base sequence, using alignment of sequences\n');
        fprintf('           ');

        for m = 1:length(j),
          fprintf(' %4s  ', t{j(m),1});
        end
        fprintf('\n');

        for m = 1:length(j),
          T = sprintf('%4s', t{j(m),1});
          T = [T sprintf('%5.2f ', n(j(m),1))];
          for nn = 1:length(j),
            if nn == m || prop(j(m),j(nn)) == 0,
              T = [T sprintf('       ')];
            else
              T = [T sprintf(' %5.1f ', 100*full(prop(j(m),j(nn))))];
            end
          end
          fprintf('%s\n',T);
        end
        fprintf('\n');

        fprintf('Sequences of longest chains\n');
        for m = 1:length(j),
          if length(t{j(m),11}) < 10000,
            fprintf('%4s %s\n', t{j(m),1}, t{j(m),11});
          end
        end
        fprintf('\n');

        fprintf('Geometric discrepancy between aligned bases\n');
        fprintf('           ');

        for m = 1:length(j),
          fprintf(' %4s  ', t{j(m),1});
        end
        fprintf('\n');

        for m = 1:length(j),
          T = sprintf('%4s', t{j(m),1});
          T = [T sprintf('%5.2f ', n(j(m),1))];
          for nn = 1:length(j),
            discr = full(Discrepancies(j(m),j(nn)));
            if nn == m || discr == Inf,
              T = [T sprintf('       ')];
            else
              T = [T sprintf(' %5.2f ', discr)];
            end
          end
          fprintf('%s\n',T);
        end
      end
    end

    ff = 1;                               % default file to use

    if Criterion > 10,                    % replace with preferred name
      pp = [];
      for f = 1:length(j),
        p = find(ismember(Preferred,t{j(f),1}));
        if ~isempty(p),
          pp = p;
          ff = f;
          fprintf('Found %s in the list of preferred structures\n', t{j(f),1});
          cn = t{j(f),1};
        end
      end
    end

    for m = 1:length(j),
      t{j(m),10} = t{j(ff),1};            % best file represents all others
    end

    c = c + 1;

    if length(j) == 1,
      Text{c} = sprintf('Unique structure is %4s, which ', t{j(ff),1});
    else
      Text{c} = sprintf('Chosen structure is %4s, which ', t{j(ff),1});
    end

    Text{c} = [Text{c} sprintf('has %4d nucleotides, %4d in longest chain, %4d pairs, ', n(j(ff),2), n(j(ff),4), n(j(ff),3))];
    Text{c} = [Text{c} sprintf('resolution %5.2f ', n(j(ff),1))];

    Text{c} = [Text{c} sprintf(' %10s | %s | %s | %s', t{j(ff),4}, t{j(ff),8}, t{j(ff),2}, t{j(ff),5})];

    fprintf('%s\n', Text{c});

  end
end

% ----------------------------------------------- List chosen files

fprintf('\nInformation about chosen files:\n');

for c = 1:length(Text),
  fprintf('%4s\n', Text{c});
end

% ----------------------------------------- Update PDBInfo with this list

if Criterion == 5 && length(t(:,1)) == OriginalNumStruct,

  save([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

  if exist('Dropboxroot') && (nargin > 1),
    save([DropboxRoot filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
  end

else

  fprintf('Note: updated PDBInfo was not saved\n');

end

fprintf('Non-redundant list has %d entries.\n', length(Text));
fprintf('Total elapsed time %8.6f minutes.\n', (cputime-tot)/60);

diary off

