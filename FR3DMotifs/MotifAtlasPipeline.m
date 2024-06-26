%==========================================================================
% The main entry point for updating the RNA 3D Motif Atlas.

% takes one argument:
% output_dir = full path to the directory where all resulting files
% should be stored.
% It is presumed that loops.txt is already written there.

% status 0 = success
% status 1 = failure

% example output_dir under MA_root: /Releases/IL_20121005_1233
%==========================================================================

function [status, err_msg] = MotifAtlasPipeline(output_dir)

    status  = 0;  % success
    err_msg = '';

    %try

        startLogging(output_dir);

        % loops.txt is already written to the output directory
        loop_id_file = fullfile(output_dir,'loops.txt')

        if ischar(loop_id_file)
            loop_ids = parse_loop_id_file(loop_id_file);
        elseif ~iscell(loop_id_file)
            status = 1;
            err_msg = 'Wrong list of loop ids';
            disp(err_msg);
            return;
        end

        disp('aCreateMM: Loading search data for all loops');
        MM = aCreateMM(loop_ids, output_dir);
        fprintf('MM has %d rows\n', size(MM,1))

        disp('aAnalyzeExtraNucleotides: Analyzing non-matched nucleotides');
        MM = aAnalyzeExtraNucleotides(MM, loop_ids, output_dir);
        fprintf('MM has %d rows\n', size(MM,1))

        disp('aSymmetrizeMatrix: Use best matching direction')
        MM = aSymmetrizeMatrix(MM, loop_ids, output_dir);
        fprintf('MM has %d rows\n', size(MM,1))

        disp('aMaximumCliques: Identify maximal cliques to form motif groups')
        groups = aMaximumCliques(MM, loop_ids, 1);

        disp('groupsToSearches: ')
        groupsToSearches(output_dir, groups);

        disp('groupsToGraphML: ')
        groupsToGraphML(output_dir, groups, MM, loop_ids, 1);

        disp('folderToVarna: ')
        folderToVarna(output_dir);

        disp('exportMotifRelease: ')
        exportMotifRelease(output_dir);


    %catch err
    %    err_msg = 'Error in MotifAtlasPipeline';
    %    disp(err_msg);
    %    status = 2;
    %end

    stopLogging();

end

function [loop_ids] = parse_loop_id_file(filename)

    loop_ids = {};

    fid = fopen(filename, 'r');

    line = fgetl(fid);
    loop_ids = regexp(line, ',', 'split');

    fclose(fid);

end

function startLogging(location)

    % log to the same directory as the output files
    % that increases the size of the output directory,
    % but keeps the log files separate

    filename = fullfile(location, 'MotifAtlasPipeline.log');

    fopen(filename, 'a');
    diary(filename);

    disp(filename);

end

function stopLogging()

    diary off;

end
