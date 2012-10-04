function [] = setup()

    % set the environmental variable with the path to the MotifAtlas folder
    motifatlas_path = [fileparts(which(mfilename)) filesep 'MotifAtlas'];
    if ~exist(motifatlas_path, 'dir'), mkdir(motifatlas_path); end        
    % set environmental variable
    setenv('MA_root', motifatlas_path);
            
    addpath(genpath('FR3DSource'));
    addpath(genpath('FR3DDevelopment'));
    addpath(genpath('FR3DMotifs'));
    addpath(genpath('PrecomputedData'));
    addpath(genpath('PDBFiles'));

end