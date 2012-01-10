function [s] = aImplode(cell_array)

    s = '';
    for i = 1:length(cell_array)
        s = [s ', ' cell_array{i}]; %#ok<AGROW>
    end
    s = s(3:end);

end
