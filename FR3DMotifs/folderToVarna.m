function [] = folderToVarna(Location)

    files = dir( fullfile(Location, '*.mat') );
    N = length(files);
    
    filename = fullfile(pwd, ['varna_' strrep(Location, filesep, '-') '.txt']);

    fid = fopen(filename,'w');

    destination = fullfile(Location, 'html');
    
    for i = 1:N
        
        clear Search;
        load( fullfile(Location, files(i).name) );
        if ~exist('Search','var')
            continue;
        end
        Search.FileName = files(i).name;
        command = motifToVarna(Search, destination, 0);
        
        fprintf(fid, '%s\n', command);
        
    end

    fclose(fid);

    unix(sprintf('bash %s',filename));
    
    movefile(filename,Location);

end