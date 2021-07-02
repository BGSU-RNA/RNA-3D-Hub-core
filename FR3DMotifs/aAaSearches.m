function [] = aAaSearches(loop_ids, start, stop, exactSizeLimit)

    if nargin < 4
        exactSizeLimit = 0;
    end

    if ischar(loop_ids) % assume file
        fid = fopen(loop_ids, 'r');
        line = fgetl(fid);
        loop_ids = regexp(line, ',','split');
        fclose(fid);
    elseif ~iscell(loop_ids)
        error('Incorrect input');
    end

    N = length(loop_ids);

    if nargin < 2
        start = 1;
        stop  = N;
    end
    if stop > N
        stop = N;
    end

    for i = start:stop

        fprintf('aAaSearches.m working on loop %4d, start=%4d, stop=%4d\n', i, start, stop);
        destination = [getSearchFolder filesep loop_ids{i}];
        if ~exist(destination,'dir'), mkdir(destination); end

        % read log file with negative results
        log = [destination filesep 'No_candidates.txt'];
        if exist(log,'file')
            done = textread(log,'%s');
        else
            done = {};
        end

        % load only once
        load(getPrecomputedDataAddress(loop_ids{i}));

        for j = 1:N
            if isempty(intersect(done, loop_ids{j}))
                pairwiseSearch(File, loop_ids{j}, exactSizeLimit);
            end
        end

    end

end
