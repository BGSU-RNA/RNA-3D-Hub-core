function [] = setup()

    current_dir = fileparts(which(mfilename));

    % set the environmental variable with the path to the MotifAtlas folder
    motifatlas_path = fullfile(current_dir, 'MotifAtlas');
    if ~exist(motifatlas_path, 'dir'), mkdir(motifatlas_path); end
    setenv('MA_root', motifatlas_path);

    fr3d_dir = fullfile(current_dir, 'FR3D');
    cd(fr3d_dir);
    addpath(current_dir);
    addpath(genpath(fullfile(current_dir, 'FR3DMotifs')));
    addpath(fullfile(fr3d_dir, 'FR3DSource'));
    addpath(fullfile(fr3d_dir, 'PrecomputedData'));
    addpath(fullfile(fr3d_dir, 'PDBFiles'));

end