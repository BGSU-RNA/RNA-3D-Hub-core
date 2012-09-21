%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input: pdb id
% output: filename with output, exit status and error message
% exit codes: 
%   0 = success
%   1 = failure
%   2 = no nucleotides in pdb file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FILENAME, status, err_msg] = getInteractions(pdb_id)

    try
        
        FILENAME = 'PairwiseInteractions.csv';
        status   = 2;
        err_msg  = '';        
                
        File = zAddNTData(pdb_id);

        if isempty(File.NT)
            return;
        end
        
        fid = fopen(FILENAME,'w');

        processMatrix(@zEdgeText, File.Edge);

        processMatrix(@zBasePhosphateText, File.BasePhosphate);

        processMatrix(@zBaseRiboseText, File.BaseRibose);    
                
        fclose(fid);
        status = 0;        
        
    catch err
        err_msg = sprintf('Error "%s" in getInteractions on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);
        status = 1;        
    end


    function processMatrix(functionHandle, matrix)

        interactions = find(matrix);

        for i = 1:length(interactions)

            [nt1, nt2] = ind2sub(File.NumNT, interactions(i));

            nt_id1 = aGetNTId(File, nt1);
            nt_id2 = aGetNTId(File, nt2);        

            textAnnotation = functionHandle( matrix(nt1, nt2) );

            if strcmp(textAnnotation,'---- ') || strcmp(textAnnotation,'----')
                continue;
            end

            fprintf(fid, '"%s","%s","%s"\n', nt_id1, nt_id2, textAnnotation);

        end

    end
    
end