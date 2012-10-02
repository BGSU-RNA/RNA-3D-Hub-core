%==========================================================================
% The main entry point for updating the RNA 3D Motif Atlas.

% takes two arguments:
% loop_ids = full path to the file with loop ids to be clustered
% location = full path to the directory where all resulting files 
% should be stored.

% status 0 = success
% status 1 = failure

% location = '/Users/anton/FR3D/MotifAtlas/Releases/ilaug'
%==========================================================================

function [status, err_msg] = MotifAtlasPipeline(loop_ids, location)

    status  = 0;  % success
    err_msg = '';
    
    try 

        startLogging();
        disp('test logging');
        return;
        if ~exist(location, 'dir'), mkdir(location); end

        if ischar(loop_ids)
            loop_ids = parse_loop_id_file(loop_ids);
        elseif ~iscell(loop_ids)           
            status = 1;
            err_msg = 'Wrong output';
            disp(err_msg);
            return;
        end
        
        MM = aCreateMM(loop_ids);
        
        MM = aAnalyzeExtraNucleotides(MM, loop_ids);

        MM = aSymmetrizeMatrix(MM, loop_ids);

        groups = aMaximumCliques(MM, loop_ids, 1);

        groupsToSearches(location, groups);

        groupsToGraphML(location, groups, MM, loop_ids, 1);    

        folderToVarna(location);

        exportMotifRelease(location);
    
        movefile([pwd filesep 'MM*.mat'], location);
        movefile([pwd filesep 'MM*.txt'], location);        
                        
    catch err
        err_msg = 'Error in createMM';
        disp(err_msg);
        status = 2;        
    end

end

function [loop_ids] = parse_loop_id_file(filename)

    loop_ids = {};
    
    fid = fopen(filename, 'r');
    
    line = fgetl(fid);    
    loop_ids = regexp(line, ',', 'split');
    
    fclose(fid);

end

function [] = startLogging()

    ma_root = getenv('MA_root');

    if ~strcmp(ma_root, '')
        log_path = fullfile(ma_root, 'logs');        
    else
        log_path = pwd;        
    end
   
    filename = fullfile(log_path, ...
                        [datestr(now, 'yyyy-mm-dd_HH-MM') ...
                        '_mlmotifatlas.txt']);

    fopen('filename', 'w');
    diary on;                    
    diary(filename);
    
    disp(filename);
    
end