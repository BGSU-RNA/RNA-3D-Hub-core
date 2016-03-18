function [status,err_msg] = aSaveLoops(Loops, save_location)

    try
        status  = 0;
        err_msg = '';

        % save loops from each pdb in a separate folder
        save_location = fullfile(save_location, Loops(1).PDBID);
        if ~exist(save_location,'dir'), mkdir(save_location); end

        N = length(Loops);

        for i = 1:N
            File = Loops(i);
            filename = fullfile(save_location, File.Filename);
            save(filename,'File');
        end

    catch err
        err_msg = sprintf('Error "%s" in aSaveLoops on line %i\n', err.message, err.stack.line);
        disp(err_msg);
        status = 1;
    end

end
