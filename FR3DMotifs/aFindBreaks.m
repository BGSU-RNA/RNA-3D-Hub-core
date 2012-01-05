function [positions] = aFindBreaks(File)

    connected = ones(1,File.NumNT);
    for i = 1:File.NumNT
        if isempty(find(File.Covalent(i,:)>0)) %#ok<EFIND>
            connected(i) = 0;
        end
    end
    
    positions = find(connected(1:end-1)==0);

end