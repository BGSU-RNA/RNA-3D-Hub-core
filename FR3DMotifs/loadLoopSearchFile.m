%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Used in conjunction with LoopSearchesLoader.py to import information
% about pairwise all-against-all searches into the database.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filename, err_msg] = loadLoopSearchFile(input_folder)
    
    filename = '';
    err_msg = '';
    
    try

        filename = fullfile(input_folder, 'SearchResults.csv');
        fid = fopen(filename, 'w');        
        
        files = dir([input_folder filesep '*.mat']);
                
        for f = 1:length(files)
            
            load(fullfile(input_folder, files(f).name));
            
            disc = Search.Discrepancy;
            if disc < 0.001
                disc = 0;
            end

            loop_id1 = Search.Query.Filename;
            loop_id2 = Search.Query.SearchFiles;
            
            % figure out which pdb file each loop came from
            if f == 1
                load(getPrecomputedDataAddress(loop_id1));
                loop_id1_pdb = File.PDBFilename;
            end                
            Search.Query.PDBFilename = loop_id1_pdb;
            load(getPrecomputedDataAddress(loop_id2));
            Search.File.PDBFilename = File.PDBFilename;            

            nt_list1 = '';                
            for i = 1:length(Search.Query.NT)
                nt_list1 = strcat(nt_list1, ',', aGetNTId(Search.Query, i));            
            end
            nt_list1 = nt_list1(2:end); % remove the last comma

            nt_list2 = '';
            for i = Search.Candidates(1,1:end-1)
                nt_list2 = strcat(nt_list2, ',', aGetNTId(Search.File, i));
            end
            nt_list2 = nt_list2(2:end); % remove the last comma
            
            fprintf(fid,'"%s","%s","%f","%s","%s"\n', loop_id1, loop_id2, disc, nt_list1, nt_list2);
            
        end
        
        fclose(fid);
                             
    catch err
        err_msg = sprintf('Error "%s" on line %i\n', err.message, err.stack.line);
        disp(err_msg);      
    end    
            
end