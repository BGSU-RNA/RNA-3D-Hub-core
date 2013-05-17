%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input: resolution ('1.5','2.0','2.5','3.0','3.5','4.0','20.0','all')
% output: filename with output, exit status and error message
% XXXX, Group1
% YYYY, Group2
% ZZZZ, Group2
%
% exit codes:
%   0 = success
%   1 = failure
%
% adapted from zWriteHTMLFielList.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filename, status, err_msg] = getNrList(resolution)

try
    status   = 0;
    err_msg  = '';
    filename = fullfile(pwd, ['nr_' resolution '.txt']);

    % set output type
    if strcmp(resolution, 'all')
        MaxRes = Inf;
    else
        MaxRes = str2double(resolution);
    end

    % load mat file with the data prepared by NR pipeline
    load PDBInfo;

    % get equivalence classes
    if MaxRes == Inf,
        % allow NMR, have at least one basepair
        i = find((n(:,2) > 0));
    else
        % Not NMR, res to MaxRes, at least one nucleotide                            
        i = find((n(:,1) > 0) .* (n(:,1) <= MaxRes) .* (n(:,2) > 0));                
    end

    [tt,nn] = zEquivalents(t(i,:),n(i,:),5,1); % get equivalents

    Keep = [];
    Keep(1) = 1;
    for j = 2:length(tt(:,1))
        if ~strcmp(tt{j,10},tt{j-1,10}) % first instance from group
            Keep(j) = 1;
        end
    end

    i = find(Keep);
    tt = tt(i,:);
    nn = nn(i,:);

    % create output files
    fid = fopen(filename, 'w');
    group = 1;
    for r = 1:length(tt(:,1))

        % representative
        fprintf(fid, '%s,Group_%i\n', tt{r,1}, group);
        
        % equivalent structures, regardless of resolution
        g = find(ismember(t(:,10),tt{r,10}));
        if length(g) > 1
            [y,i] = sort(-n(g,3)./n(g,2)); % sort by pairs per nucleotide
            g = g(i);
            for h = 1:length(g)
                if ~strcmpi(t{g(h),1},tt{r,1})
                    fprintf(fid, '%s,Group_%i\n', t{g(h),1}, group);
                end
            end
        end            
        group = group + 1;
        
    end
    fclose(fid);

catch err
    err_msg = sprintf('Error "%s" in getNrList on line %i\n', err.message, err.stack.line);
    disp(err_msg);
    status = 1;
end

end