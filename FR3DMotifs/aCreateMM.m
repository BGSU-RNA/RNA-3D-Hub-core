function [MM] = aCreateMM(loop_ids,output_dir)

    tic;

    FILENAME = fullfile(output_dir,'MM_initial.mat');

    % if the file has already been created, use that
    if exist(FILENAME,'file')
        disp('Loading MM matrix from MM_initial.mat')
        load(FILENAME,'MM');
        return
    end

    N = length(loop_ids);

    SavedSearchLocation = getSearchFolder;

    NO_CANDIDATES_FILE = 'No_candidates.txt';

    NO_MATCH = 10;
    INITIAL_VALUE = -1;

    MM(1:N,1:N) = INITIAL_VALUE;


    for i = 1:N

        fprintf('aCreateMM: Loading search data for %s, loop %i out of %i\n', loop_ids{i}, i, N);

        subdir = fullfile(SavedSearchLocation, loop_ids{i});

        process_no_candidates_file;
        load_search_data;
        % save_intermediate_results;

    end


    save(FILENAME, 'MM', 'loop_ids');

    check_never_updated;

    % diagnostic information, printed but not acted on
    checkMatchingMatrix(MM);

    toc;


    % NO_CANDIDATES_FILE stores a list of loop_ids that didn't match
    % this loop_id during all-against-all searches. This function reads
    % the file and sets the corresponding cells to NO_MATCH
    function process_no_candidates_file

        logfile = fullfile(subdir, NO_CANDIDATES_FILE);

        if exist(logfile,'file')

        	fid = fopen(logfile, 'r');
        	no_candidates = textscan(fid, '%s');
        	fclose(fid);
            [a,b] = intersect(loop_ids, no_candidates{1});
            b = reshape(b, 1, []);
            if ~isempty(b)
                MM(i,b) = NO_MATCH;
            end

        end

    end


    % List all mat files in the directory, identify which files need to be
    % loaded, load them and store the discrepancy in the matrix.
    function load_search_data

        filelist = dir( fullfile(subdir, '*.mat') );

        if isempty(filelist), return; end

        targets = arrayfun(@(x) x.name(13:23), filelist, 'UniformOutput', false);

        [a,b] = intersect(loop_ids, targets);
        b = reshape(b, 1, []);

        for j = 1:length(b)

            search_file = fullfile(subdir, strcat(loop_ids{i}, '_', loop_ids{b(j)}, '.mat'));
            try
                load(search_file, 'Search');
            catch
                fprintf('Corrupted file %s\n', search_file);
                pairwiseSearch(loop_ids{i}, loop_ids{b(j)});
                load(search_file, 'Search');
            end

            disc = min(Search.Discrepancy);
            if isempty(disc)
                MM(i,b(j)) = NO_MATCH;
            else
                % set discrepancy to 0 if name is the same
                if i == b(j) && disc < 0.1
                    MM(i,b(j)) = 0;
                else
                    MM(i,b(j)) = disc;
                end
            end

        end
    end


    % The matrix is initialized with INITIAL_VALUE. This value should be
    % updated either with NO_MATCH, or with Discrepancy, or with 0.
    function check_never_updated

        [r, c] = find( MM(i,:) == INITIAL_VALUE );
        for j = 1:length(r)
            MM(r(j), c(j)) = pairwiseSearch( loop_ids{r(j)}, loop_ids{c(j)} );
        end

        if ~isempty(r)
            save(FILENAME,'MM','loop_ids');
        end

    end


    % save every 50 iterations
    function save_intermediate_results
        if ~mod(i, 50)
            save(FILENAME,'MM','loop_ids');
        end
    end

end
