function [] = aAaSearches(loop_ids, start, stop)


    N = length(loop_ids);
    
    SearchResults = getSearchFolder;

    if nargin < 2
        start = 1;
        stop  = N;
    end
    if stop > N
        stop = N;
    end
    
    for i = start:stop    
        
        disp(i);
        destination = [SearchResults filesep loop_ids{i}];
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
            % don't repeat the search if 
            if isempty(intersect(done, loop_ids{j}))            
                pairwiseSearch(File, loop_ids{j});
            end
        end
        
    end
    
end