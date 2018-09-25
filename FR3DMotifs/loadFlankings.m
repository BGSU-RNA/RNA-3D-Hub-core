%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input: pdb id
% output: filename with output, exit status and error message
% exit codes:
%   0 = success
%   1 = failure
%   2 = no nucleotides in pdb file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FILENAME, status, err_msg] = loadFlankings(pdb_id)

    try

        FILENAME = fullfile(pwd, 'FlankingInteractions.csv');
        status   = 2;
        err_msg  = '';

        File = zAddNTData(pdb_id);

        if isempty(File.NT)
            return;
        end

        fid = fopen(FILENAME,'w');

        processMatrix(File.Flank);

        fclose(fid);
        status = 0;

    catch err
        err_msg = sprintf('Error "%s" in loadFlankings on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);
        status = 1;
    end


    function processMatrix(matrix)

        interactions = find(matrix);

        for i = 1:length(interactions)

            [nt1, nt2] = ind2sub(File.NumNT, interactions(i));

            nt_id1 = File.NT(nt1).ID;
            nt_id2 = File.NT(nt2).ID;
           

            fprintf(fid, '"%s","%s","%i"\n', nt_id1, nt_id2, 1);

        end

    end

end
