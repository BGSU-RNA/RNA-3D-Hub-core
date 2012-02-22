function search_results_equivalence_test()

    folder1 = 'MotifAtlas/aAa';
    folder2 = 'SearchDB';
    
    list = dir(folder1);
    c= 0;
    d = 0;
    for i = 1:length(list)
        if list(i).name(1) == '.'
            continue;
        end
        list2 = dir(fullfile(folder1, list(i).name, ''));
        for j = 1:length(list2)
            if ~exist(fullfile(folder2, list(i).name, list2(j).name), 'file')
                disp(fullfile(folder2, list(i).name, list2(j).name));
                c = c + 1;
            else
                fprintf('found in both\n');
                d = d + 1;
            end            
        end
        
        % test the reverse
        list3 = dir(fullfile(folder2, list(i).name,''));
        for j = 1:length(list3)
            if list3(j).name(1) ~= '.' && ...
               ~strcmp(list3(j).name,'No_candidates.txt') && ...
               strcmp(list3(j).name(16:19),'1J5E')
                if ~exist(fullfile(folder1, list(i).name, list3(j).name), 'file')
                    disp(fullfile(folder1, list(i).name, list3(j).name))
                    keyboard;
                else
                    fprintf('found in both\n');
                end
            end
        end 

    end
    c
    d

    
end