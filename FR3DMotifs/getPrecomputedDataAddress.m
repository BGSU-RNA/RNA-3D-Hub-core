function [address] = getPrecomputedDataAddress(motif)
    
    address = fullfile(getenv('MA_root'), 'PrecomputedData', motif(4:7), [motif '.mat'] );

end