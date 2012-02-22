% pairwise search address, loop - loop search

function [address] = getSearchAddress(motif1, motif2)
        
    address = fullfile(getSearchFolder, motif1, [motif1 '_' motif2 '.mat']);
       
end