function [Search] = aFR3DSearch(Query,GenericFile)

Search = struct;
Filenames = Query.SearchFiles;
Verbose = 1;

if nargin < 2
    [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
else
    File = GenericFile.File;
    SIndex = GenericFile.QIndex;
end    
% ------------------------------------------- Store actual filenames
%                                             rather than list name(s)

clear Filenames;

for i = 1:length(SIndex),
    Filenames{i} = File(SIndex(i)).Filename;
end

% ------------------------------------------- Construct details of search
if isfield(Query,'Filename'),                % if query motif is from a file
    [File,QIndex] = zAddNTData(Query.Filename,0,File);  
    Query = xConstructQuery(Query,File(QIndex)); % preliminary calculations
else
    Query = xConstructQuery(Query);              % preliminary calculations
end


disp(Query);
clear Search
Search.SaveName = [datestr(now,31) '-' Query.Name];  
                                  % use date and time to identify this search

% ------------------------------------------- Display query information------


  fprintf('Query %s:', Query.Name);             % display query name

  if isfield(Query,'Description'),
    fprintf(' %s\n', Query.Description);
  else
    fprintf('\n');
  end
% ------------------------------------------- Calc more distances if needed -

for f=1:length(SIndex),
  i = SIndex(f);
  if isempty(File(i).Distance),
    dmin = 0;
  else
    dmin = ceil(max(max(File(i).Distance)));
  end

  if (ceil(Query.DistCutoff) > dmin) && (File(i).NumNT > 0),
    c = cat(1,File(i).NT(1:File(i).NumNT).Center);
    File(i).Distance = zMutualDistance(c,Query.DistCutoff); 
  end
end

drawnow

% ------------------------------------------- Find candidates ---------------

starttime = cputime;
Candidates = xFindCandidates(File(SIndex),Query,Verbose);  % screen for candidates

if ~isempty(Candidates),                         % some candidate(s) found
    
     if Query.Geometric > 0,
          [Discrepancy, Candidates] = xRankCandidates(File(SIndex),Query,Candidates);
          fprintf('Found %d candidates in the desired discrepancy range\n',length(Discrepancy));

           if (Query.ExcludeOverlap > 0) && (~isempty(Discrepancy)) ...
             && (Query.NumNT > 2),
             [Candidates, Discrepancy] = xReduceOverlap(Candidates,Discrepancy); 
                                                         % quick reduction in number
             [Candidates, Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,400); 
                                                        % find top 400 distinct ones
             fprintf('Removed highly overlapping candidates, kept %d\n', length(Candidates(:,1)));
           end

     elseif Query.NumNT > 2,
          A = [Candidates sum(Candidates')'];        % compute sum of indices
          N = Query.NumNT;                           % number of nucleotides
          [y,i] = sortrows(A,[N+1 N+2 1:N]);         % sort by file, then this sum
          Candidates = Candidates(i,:);              % put all permutations together
          Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates
     else
          N = Query.NumNT;                           % number of nucleotides
          [y,i] = sortrows(Candidates,[N+1 1 2]);
          Candidates = Candidates(i,:);              % put all permutations together
          Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates
     end

    % -------------------------------------------------- Save results of search

      Search.Query       = Query;
      Search.Filenames   = Filenames;
      Search.TotalTime   = cputime - starttime;
      Search.Date        = Search.SaveName(1:10);
      Search.Time        = Search.SaveName(12:18);
      Search.SaveName    = strrep(Search.SaveName,' ','_');
      Search.SaveName    = strrep(Search.SaveName,':','_');
      Search.SaveName    = strrep(Search.SaveName,'<','_');
      Search.SaveName    = strrep(Search.SaveName,'>','_');
      Search.SaveName    = strrep(Search.SaveName,'?','_');
      Search.SaveName    = strrep(Search.SaveName,'*','_');
      Search.SaveName    = strrep(Search.SaveName,'&','_');
      Search.Candidates  = Candidates;
      Search.Discrepancy = Discrepancy;

      Search = xAddFiletoSearch(File(SIndex),Search);

    % ------------------------------------------------ Display results
    fprintf('Entire search took %8.4f seconds, or %8.4f minutes\n', (cputime-starttime), (cputime-starttime)/60);

    if ~isempty(Candidates)
        xListCandidates(Search,Inf,1);
        if isfield(Query,'SaveDir') 
           outdir = Query.SaveDir;
        else
           outdir = 'tempResults'; 
        end

        if (~exist(outdir,'dir')),
            mkdir(outdir);
        end     
        save([outdir filesep Query.Name], 'Search');
    end

end

end