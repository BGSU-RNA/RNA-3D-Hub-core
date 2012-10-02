%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Used in conjunction with LoopSearchesLoader.py to import information
% about nucleotide positions in the loop into the database.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filename, err_msg] = loadLoopPositions(input_folder)
    
    filename = '';
    err_msg = '';
    
    try

        filename = fullfile(input_folder, 'LoopPositions.csv');
        fid = fopen(filename, 'w');        
        
        files = dir([input_folder filesep '*.mat']);
                
        for f = 1:length(files)
            
            load(fullfile(input_folder, files(f).name));
            
            loop_id = File.Filename;
            bulges = aDetectBulgedBases(File);
            
            for i = 1:length(File.NT)
                position = i;
                nt_id = aGetNTId(File, i);
                if ~isempty(find(full(File.Flank(i,:)),1)),
                    flank = 1;
                else
                    flank = 0;
                end
                if ~isempty(find(bulges == i, 1))
                    bulge = 1;
                else
                    bulge = 0;
                end                                
                fprintf(fid,'"%s","%i","%s","%i","%i"\n', loop_id, position, nt_id, bulge, flank);
            end
        end
        
        fclose(fid);
                             
    catch err
        err_msg = sprintf('Error "%s" on line %i\n', err.message, err.stack.line);
        disp(err_msg);      
    end    
            
end