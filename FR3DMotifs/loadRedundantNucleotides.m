%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Used in conjunction with RedundantNucleotidesLoader.py to import information
% about redundant nucleotides within pdb files into the database.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filename, err_msg] = loadRedundantNucleotides(pdb_id)
    
    filename = '';
    err_msg = '';
    
    try

        filename = fullfile(pwd, 'RedundantNucleotides.csv');
        fid = fopen(filename, 'w');        
        
        F = zAddNTData(pdb_id);
                
        if ~isfield(F, 'Redundant'), return; end        
                
        indexes = find(F.Redundant);
        N = length(indexes);
        for i = 1:N
            [x, y] = ind2sub(F.NumNT, indexes(i));
            nt_id1 = aGetNTId(F, x);
            nt_id2 = aGetNTId(F, y);
            fprintf(fid, '"%s","%s"\n', nt_id1, nt_id2);            
        end               
        
        fclose(fid);
                             
    catch err
        err_msg = sprintf('Error "%s" on line %i\n', err.message, err.stack.line);
        disp(err_msg);      
    end    
            
end