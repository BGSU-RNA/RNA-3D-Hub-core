
load PDBInfo

reportdate = '2011-05-21';

blackandwhite = 1;                       % bw graphs and fewer of them, if 1

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

Timeline = [];

if ~exist('reportdate')
  reportdate = datestr(date,'yyyy-mm-dd');
end

tot = cputime;                    % keep track of total time to run

Criterion = 5;                   
                                  % 1-earliest release date 
                                  % 2-resolution 
                                  % 3-number of nucleotides
                                  % 4-#pairs
                                  % 5-highest #BP / #nucleotides
                                  % 6-most recent release date

                                  % add 10, use preferred list to override

% --------------------------------- Look up information on these structures

fprintf('Starting with %4d structures from the PDB.\n', length(t(:,1)));

for i = 1:length(t(:,1)),
  n(i,4) = length(t{i,11});       % length of longest chain
end

if 0 > 1,
  i = find(n(:,3) > 0);             % restrict to structures with a basepair
  t = t(i,:);
  n = n(i,:);           
  fprintf('Found %4d structures with at least one basepair.\n', length(t(:,1)));
end

if 10 > 1,
  i = find(n(:,2) > 0);             % restrict to structures with a nucleotide
  t = t(i,:);
  n = n(i,:);           
  fprintf('Found %4d structures with at least one nucleotide.\n', length(t(:,1)));
end

[y,i] = sort(n(:,4));         % sort by number of nucleotides in longest chain
t = t(i,:);
n = n(i,:);

F = length(i);                    % number of files

clear d

for i = 1:F,
  depdate = t{i,4};
  if length(depdate) < 8,
    depdate = ['0' depdate];
  end

  d(i,1) = 1995 + (datenum(depdate,'mm-dd-yy') - datenum('01-01-95','mm-dd-yy'))/365;  % dep date

  if d(i,1) < 1000,
    depdate
  end

end

now = 1995 + (datenum('03-26-11','mm-dd-yy') - datenum('01-01-95','mm-dd-yy'))/365; 

figure(4)
clf

if blackandwhite > 0,
  ResolutionLimits = [Inf 1 2 3 4 100];   % Inf must be first!
else
  ResolutionLimits = [Inf 1 1.5 2 2.5 3 3.5 4 100];   % Inf must be first!
end

for r = 1:length(ResolutionLimits),      % loop through resolution limits

MaxRes = ResolutionLimits(r);

if r == 2,
  n = nn;                                    % store from first pass
end

if MaxRes == Inf,
  tt = t;
  nn = n;
  dd = d;
else
  i = find(n(:,1) <= MaxRes);
  tt = t(i,:);
  nn = n(i,:);
  dd = d(i,:);
end

F = length(dd);                              % number of structures left

fprintf('\n');
fprintf('Found %d structures with resolution at or better than %7.4f #################\n',F,MaxRes)

% -------------------------------------------

fprintf('Structures are sorted by ');
switch mod(Criterion,10)
  case 1, fprintf('release date, earliest first\n');
  case 2, fprintf('resolution, best resolution first\n');
  case 3, fprintf('number of nucleotides (decreasing)\n');
  case 4, fprintf('number of basepairs (decreasing)\n');
  case 5, fprintf('number of basepairs per nucleotide (decreasing)\n');
  case 6, fprintf('release date, most recent first\n');
end

done = zeros(1,F);                % whether each file has been considered
c = 0;                            % counts the number of files selected so far

Timeline = [];                    % start with an empty timeline

if F == 0,
  Timeline = [[now-10 0 0 0]; [now 0 0 0]];
end

for i = 1:(F-1),                  % loop through all files at this resolution
  if done(i) == 0,                % file does not already appear in a report

    j = find(ismember(tt(:,10),tt{i,10})); % equivalence class i lies in

    if length(j) == 1,             % no other file is equivalent to i
      done(i) = 1;                 % no need to display this one again
      crit(1,Criterion) = 0;       % good value of criterion, for later
      y(1) = dd(j(1),1);           % deposition date
    elseif length(j) == 0,
      disp('Trouble!');            % do nothing, but this case shouldn't happen
    else
      crit = [];                   % criterion for sorting and choosing

      for f = 1:length(j),         % go through the equivalence class
        done(j(f)) = 1;            % no need to display this one again
                                   % sorting criteria
        crit(f,1) = dd(j(f),1);           % deposition date
        crit(f,2) = nn(j(f),1);            % resolution
        crit(f,3) = nn(j(f),2);           % number of nucleotides
        crit(f,4) = nn(j(f),3);               % number of pairs
        crit(f,5) = nn(j(f),3)/nn(j(f),2); % ratio # bp / # nucleotides
      end

%      fprintf('\n');

       [y,k] = sortrows(crit,[1 -5 2 3 4]);% sort equiv class by deposition date

      j = j(k);                         % reorder class by deposition date
      crit = crit(k,:);                 % reorder criterion data too!
    end

    ff = 1;                               % default file to use

    maxcrit = -1;                         % current best value of criterion
    curNT = 0;                            % current number of nucleotides
    curNP = 0;                            % current number of pairs
    curFile = '    ';                     % current best file

    for m = 1:length(j),                  % go through the equivalence class
      if crit(m,Criterion) > maxcrit,     % a better file was found

        newFile = tt{j(m),1};              % store file name

fprintf('Group %4s file %4s resolution %4.1f was added on %10s (%7.6f) with %4d nucleotides, %4d basepairs, ratio %7.4f\n', tt{j(m),10}, tt{j(m),1}, nn(j(m),1), tt{j(m),4}, y(m), nn(j(m),7), nn(j(m),8), nn(j(m),3)/nn(j(m),2));

        newNT = nn(j(m),7);
        newNP = nn(j(m),8);

        if m == 1,
          Timeline = [Timeline; [y(m) 1 (newNT-curNT) (newNP-curNP)]]; % new class
        else
          Timeline = [Timeline; [y(m) 0 (newNT-curNT) (newNP-curNP)]]; % new file in class
        end

        if newNT < curNT || newNP < curNP,
          fprintf('Reduction of %4d nucleotides %4d pairs between %4s and %4s\n', newNT-curNT, newNP-curNP, curFile, newFile);
        end

        curNT = newNT;
        curNP = newNP;
        curFile = newFile;
        maxcrit = crit(m,Criterion);
      end
    end
    Timeline = [Timeline; [now 0 0 0]];
  end
end

% ---------------------------------------------- Graphs of growth of database

  [y,i] = sort(Timeline(:,1));                 % sort by date of increase
  Timeline = Timeline(i,:);                    % re-order data

  Date  = Timeline(:,1);
  NumF  = Timeline(:,2);
  NumNT = Timeline(:,3);
  NumP  = Timeline(:,4);

  Year = 1995 + (Date - datenum('01/01/1995','mm/dd/yyyy'))/365;
  Year = Date;

  color = 'krgcmbkrgcmb';

colors = [[204 0 255]/255; ...            % All  purple
          [0 0 0];     ...                % 1A   black
          0.4*[1 1 1]; ...                % 1.5A gray
          [0 0 1]; ...                    % 2A   blue
          [0 1 1]; ...                    % 2.5A cyan
          [0 1 0]; ...                    % 3A   green
          [245 245 15]/255; ...           % 3.5A yellow
          [255 165 0]/255; ...            % 4A   orange
          [1 0 0]; ...                    % 100A red
         ];

  if blackandwhite > 0,
    colors = 0*colors;                      % black and white
  end

fs = 10;

  subplot(1,2,1)
  stairs(Year,cumsum(NumF),'color',colors(r,:));
  title('Number of equivalence classes','fontsize',fs);
  hold on

  subplot(1,2,2)
  stairs(Year,cumsum(NumNT),'color',colors(r,:));
  title('Nucleotides in representative structures','fontsize',fs);
  hold on

if 0 > 1,
  subplot(1,3,3)
  stairs(Year,cumsum(NumP),'color',colors(r,:));
  title('Number of basepairs in distinct structures');
  hold on
end

end

fprintf('Total elapsed time %8.6f minutes.\n', (cputime-tot)/60);

if 10 > 1,
  subplot(1,2,1)
  axis([1993 now+1/12 0 1000]);
  hold on
  set(gca,'fontsize',fs);

  subplot(1,2,2)
  axis([1993 now+1/12 0 92000]);
  set(gca,'YTick',[0 10000 20000 30000 40000 50000 60000 70000 80000 90000])
  set(gca,'YTickLabel',{'0','10000','20000','30000','40000','50000','60000','70000','80000','90000'},'fontsize',fs);
  hold on
  set(gca,'fontsize',fs);

  if 0 > 1,
  subplot(1,3,3)
  axis([1993 now+1/12 0 35000]);
  set(gca,'YTick',[0 5000 10000 15000 20000 25000 30000 35000])
  set(gca,'YTickLabel',{'0','5000','10000','15000','20000','25000','30000','35000'},'fontsize',fs);
  hold on
  end

  orient landscape
  saveas(gcf,'growth_of_NR_set.eps');
  saveas(gcf,'growth_of_NR_set.png');
  saveas(gcf,'growth_of_NR_set.pdf');

end
