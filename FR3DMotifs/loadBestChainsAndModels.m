%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Used in conjunction with .py to import information
% about best chains and models within pdb files into the database.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [best_chains, best_models, err_msg] = loadBestChainsAndModels(pdb_id)

    best_chains = '';
    best_models = '';
    err_msg = '';

    try
        F = zAddNTData(pdb_id);

        % cell array
        best_chains = zBestChains(F);
        best_chains = strjoin(reshape(best_chains{1}, 1, []),',');

        % cell array with integer arrays
        best_models = num2str(cell2mat(zBestModels(F)), '%i,');
        best_models = best_models(1:end-1);

    catch err
        err_msg = sprintf('Error "%s" on line %i\n', err.message, err.stack.line);
        disp(err_msg);
    end

end
