function [pdb] = aParsePdbFilename(filename)

    r = regexp(filename,'[a-zA-Z0-9]{4}','match');
    if isempty(r)
        pdb = '';
    else
        pdb = r{1};
    end

end