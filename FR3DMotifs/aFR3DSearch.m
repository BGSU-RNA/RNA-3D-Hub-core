function [Search] = aFR3DSearch(Query,GenericFile)

Filenames = Query.SearchFiles;
Verbose = 0;
MaxNumberOfCandidates = Inf;
% MaxNumberOfCandidates = 1000;

if nargin < 2

    [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data

    % File.NumNT = -1 when the PDB file could not be found locally or in PDB.
    if File.NumNT < 0
        error('PDB file %s was not properly read', File.Filename);
    end

    % don't search if there are no nucleotides in the file, otherwise will crash.
    if isempty(File(1).NT)
        Search = [];
        return;
    end
    Query = xConstructQuery(Query, File);

else

    File = GenericFile.File;
    % File(1) - Query File, File(2) - target file
    Query = xConstructQuery(Query,File(1));
    File(1) = [];
    SIndex = 1;

end

Query.Verbose = 0;
fprintf('Query %s\n', Query.Name);

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

% ------------------------------------------- Find candidates ---------------

%fprintf('DEBUG(fs): before xFindCandidates\n');
starttime = cputime;
Candidates = xFindCandidates(File(SIndex),Query,Verbose);  % screen for candidates
%fprintf('DEBUG(fs): after xFindCandidates\n');


if ~isempty(Candidates)                         % some candidate(s) found

    if Query.Geometric > 0

        [Discrepancy, Candidates] = xRankCandidates(File(SIndex),Query,Candidates);
        fprintf('Found %d candidates in the desired discrepancy range\n',length(Discrepancy));

        if (Query.ExcludeOverlap > 0) && (~isempty(Discrepancy)) && (Query.NumNT >= 2)

            [Candidates, Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,MaxNumberOfCandidates);
            [Candidates, Discrepancy] = xExcludeRedundantCandidates(File(SIndex),Candidates,Discrepancy);

        end

    elseif Query.NumNT > 2,

        if (Query.ExcludeOverlap > 0)
            [Candidates] = xExcludeRedundantCandidates(File(SIndex),Candidates);
            if Verbose > 0,
                fprintf('Removed candidates from redundant chains, kept %d\n', length(Candidates(:,1)));
            end
        end

        A = [Candidates sum(Candidates')'];        %#ok<UDIM> % compute sum of indices
        N = Query.NumNT;                           % number of nucleotides
        [y,i] = sortrows(A,[N+1 N+2 1:N]);         % sort by file, then this sum
        Candidates = Candidates(i,:);              % put all permutations together
        Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates

    else

        if (Query.ExcludeOverlap > 0)
            [Candidates] = xExcludeRedundantCandidates(File(SIndex),Candidates);
            if Verbose > 0,
                fprintf('Removed candidates from redundant chains, kept %d\n', length(Candidates(:,1)));
            end
        end

        N = Query.NumNT;                           % number of nucleotides
        [y,i] = sortrows(Candidates,[N+1 1 2]);
        Candidates = Candidates(i,:);              % put all permutations together
        Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates

    end

    % -------------------------------------------------- Save results of search

    Search.SaveName    = [datestr(now,31) '-' Query.Name];
    Search.Query       = Query;
    Search.Filenames   = Filenames;
    Search.TotalTime   = cputime - starttime;
    Search.Date        = Search.SaveName(1:10);
    Search.Time        = Search.SaveName(12:18);
    Search.Candidates  = Candidates;
    Search.Discrepancy = Discrepancy;

    fprintf('DEBUG(fs): before xAddFiletoSearch\n');
    Search = xAddFiletoSearch(File(SIndex),Search);
    fprintf('DEBUG(fs): after xAddFiletoSearch\n');

    if ~isempty(Search.Candidates)
        if isfield(Query,'SaveDir')
            outdir = Query.SaveDir;
        else
            outdir = 'tempResults';
        end

        fprintf('DEBUG(fs): outdir = %s\n', outdir);

        if ~exist(outdir,'dir')
            fprinf('DEBUG(fs): creating outdir\n');
            mkdir(outdir);
        end

        fprintf('DEBUG(fs): before save attempt\n');
        Query
        Search
        [outdir filesep Query.Name]
        exist([outdir filesep Query.Name])
        save([outdir filesep Query.Name], 'Search');
        fprintf('DEBUG(fs): after save attempt\n');
    end

    fprintf('DEBUG(fs): HERE\n');

else

    fprintf('DEBUG(fs): else clause\n');
    Query
    Search.Query = Query;

end

fprintf('DEBUG(fs): end of function\n');

end
