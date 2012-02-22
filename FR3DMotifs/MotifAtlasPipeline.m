%
% 0 = success
% 1 = aAaSearches crashed
%

function [status, err_msg] = MotifAtlasPipeline(loop_ids)

    status  = 0;  % success
    err_msg = '';

    setenv('MA_root', fullfile(pwd, 'MotifAtlas', ''));
    
    try 
        aAaSearches(loop_ids);
    catch err
        err_msg = 'Error in aAaSearches';
        disp(err_msg);
        status = 1;        
    end

    try 
        createMM(loop_ids);
    catch err
        err_msg = 'Error in createMM';
        disp(err_msg);
        status = 2;        
    end








end