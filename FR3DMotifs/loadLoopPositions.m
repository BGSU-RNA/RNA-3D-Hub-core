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
            bulges  = aDetectBulgedBases(File);
            borders = [1 File.chain_breaks (File.chain_breaks + 1) File.NumNT];
            
            N = length(File.NT);
            
            for i = 1:N

                nt_id = File.NT(i).ID;
                
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

                if ~isempty(find(borders == i, 1))
                    border = 1;
                else
                    border = 0;
                end
                
                fprintf(fid,'"%s","%i","%s","%i","%i","%i"\n', loop_id, i, nt_id, bulge, flank, border);
            end
        end
        
        fclose(fid);
                             
    catch err
        err_msg = sprintf('Error "%s" on line %i\n', err.message, err.stack.line);
        disp(err_msg);      
    end    
            
end