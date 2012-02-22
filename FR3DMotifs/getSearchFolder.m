function [folder] = getSearchFolder

    folder = fullfile(getenv('MA_root'), 'aAa', '');
    
    if ~exist(folder,'dir'), mkdir(folder); end

end