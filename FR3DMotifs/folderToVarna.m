function [] = folderToVarna(location)

    groups_location = [location filesep 'Groups'];
    files = dir( fullfile(groups_location, '*.mat') );
    N = length(files);
    
    filename = fullfile(pwd, ['varna_' strrep(location, filesep, '-') '.txt']);

    fid = fopen(filename,'w');
    
    for i = 1:N
        
        clear Search;
        load( fullfile(groups_location, files(i).name) );
        if ~exist('Search','var')
            continue;
        end
        Search.FileName = files(i).name;
        command = motifToVarna(Search, location, 0);
        
        fprintf(fid, '%s\n', command);
        
    end

    fclose(fid);

    unix(sprintf('bash %s',filename));
    
    movefile(filename, location);

end