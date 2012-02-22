function [MM, names] = aCreateMM(Location, names)

    disp('Loading data...');
    
    tic;
    
%     Searches = aGetSearchFolder(Location);
    sysfilesep = filesep;
    noMatch = 10; 
    
%     names = aCreateSearchFileList(files_to_analyze,Location,1);    
    N = length(names);

    MMFileName = 'MM_initial.mat';
    
    MM(1:N,1:N) = -1;

    for i = 1:N

        fprintf('%i out of %i\n',i,N);
        subdir = [Location sysfilesep names{i}];
        
        % figure out what searches didn't match
        logfile = [subdir sysfilesep 'No_candidates.txt'];        
        if exist(logfile,'file')
            
        	fid = fopen(logfile,'r');
        	no_candidates = textscan(fid,'%s');
        	fclose(fid);
            [a,b] = intersect(names,no_candidates{1});
            if ~isempty(b)
                MM(i,b) = noMatch;
            end
            
        end
        
        % get discrepancies from matched searches
        filelist = dir([subdir sysfilesep '*.mat']);        
        if isempty(filelist)
            continue;
        end
        targets = cell(1,length(filelist));
        for j = 1:length(filelist)
            targets{j} = filelist(j).name(13:23);            
        end        
        [a,b] = intersect(names,targets);

        for j = 1:length(b)
            search_file = [subdir sysfilesep names{i} '_' names{b(j)} '.mat'];            
            try
                load(search_file);
            catch
                fprintf('Corrupted file %s\n',search_file);
                keyboard;
%                 aPairwiseSearch(names{i},names{b(j)},Location);
%                 load(search_file);                

            end
            
            Disc = min(Search.Discrepancy);
            if isempty(Disc)
                MM(i,b(j)) = noMatch;
            else
                if i == b(j) && Disc < 0.1
                    MM(i,b(j)) = 0;
                else
                    MM(i,b(j)) = Disc;
                end
            end
        end
        
        if mod(i,50) == 0 % save every 50 iterations
            save(MMFileName,'MM','names');        
        end

    end
    
%     save(MMFileName,'MM','names');        

    to_delete = [];
    for i = 1:N;
        if all(MM(i,:)<0)
            to_delete = [to_delete i]; %#ok<AGROW>
        end
    end
    
%     for i = 1:N
%         for j = 1:N
%             if MM(i,j) < 0
%                 disc = aPairwiseSearch(names{i},names{j},Location);
%                 if disc < Inf
%                     MM(i,j) = disc;
%                 end
%             end
%         end
%     end
    
    
    if ~isempty(to_delete)
    
        MM(to_delete,:) = [];
        MM(:,to_delete) = [];

        fprintf('\tNo search data about %i files (probably skipped bulges):\n',length(to_delete));
        disp(names(to_delete));
        names(to_delete) = [];
        
    end
    
    x = find(diag(MM)<0);
    MM(x,:)=[];
    MM(:,x)=[];
    names(x)=[];    
    MM(MM<0) = noMatch;

    save(MMFileName,'MM','names');            
    
    [a,b] = size(MM);
    fprintf('\tMatrix dimensions: %i by %i\n',a,b);
%     aCheckMatrix(MM);
    
    toc;

end
