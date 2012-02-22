function [bulged] = aDetectBulgedBases(File)

    bulged = [];
    for i = 1:length(File.NT)
        if isempty(find(File.Edge(i,:))) %#ok<EFIND>
            bulged(end+1) = i;
        end
    end
    
end