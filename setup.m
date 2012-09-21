function [] = setup()

    setenv('MA_root', '/Users/anton/FR3D/MotifAtlas');

    addpath(genpath('FR3DSource'));
    addpath(genpath('FR3DDevelopment'));
    addpath(genpath('FR3DMotifs'));
    addpath(genpath('PrecomputedData'));
    addpath(genpath('PDBFiles'));

end