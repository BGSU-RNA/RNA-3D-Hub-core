function MM = aCreateMM(loop_ids)

    disp('Loading data...');
        
    tic;

    N = length(loop_ids);
    
    Location = getSearchFolder;
    
    FILENAME = 'MM_initial.mat';
    NO_CANDIDATES_FILE = 'No_candidates.txt';
    
    NO_MATCH = 10;   
    INITIAL_VALUE = -1;
    
    MM(1:N,1:N) = INITIAL_VALUE;

    
    for i = 1:N

        fprintf('%i out of %i\n', i, N);
        
        subdir = fullfile(Location, loop_ids{i});
        
        process_no_candidates_file;        
        load_search_data;        
        save_intermediate_results;
        
    end
    
    
    save(FILENAME, 'MM', 'loop_ids');
    
    check_never_updated;
    
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



%         targets = cell(1, length(filelist));
%         for j = 1:length(filelist)
%             targets{j} = filelist(j).name(13:23);
%         end

%     to_delete = [];
%     for i = 1:N;
%         if all(MM(i,:)<0)
%             to_delete = [to_delete i]; %#ok<AGROW>
%         end
%     end
    
%     for i = 1:N
%         for j = 1:N
%             if MM(i,j) < 0
%                 disc = aPairwiseSearch(loop_ids{i},loop_ids{j},Location);
%                 if disc < Inf
%                     MM(i,j) = disc;
%                 end
%             end
%         end
%     end
    
    
%     if ~isempty(to_delete)
%     
%         MM(to_delete,:) = [];
%         MM(:,to_delete) = [];
% 
%         fprintf('\tNo search data about %i files (probably skipped bulges):\n',length(to_delete));
%         disp(loop_ids(to_delete));
%         loop_ids(to_delete) = [];
%         
%     end
    
%     x = find(diag(MM)<0);
%     MM(x,:)=[];
%     MM(:,x)=[];
%     loop_ids(x)=[];    
%     MM(MM<0) = NO_MATCH;

%     save(FILENAME,'MM','loop_ids');            
    
%     [a,b] = size(MM);
%     fprintf('\tMatrix dimensions: %i by %i\n',a,b);
