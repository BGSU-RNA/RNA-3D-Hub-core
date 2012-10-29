%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input: pdb id
% output: filename with output, exit status and error message
% exit codes:
%   0 = success
%   1 = failure
%   2 = no nucleotides in pdb file
% file format: "id1","id2","dist"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FILENAME, status, err_msg] = aGetDistances(pdb_id)

    try

        MAXDISTANCE = 16; % neighborhood in Angstroms
        FILENAME    = fullfile(pwd, 'Distances.csv');
        status      = 2;  % no nucleotides in pdb file
        err_msg     = '';

        F = zAddNTData(pdb_id);

        if isempty(F.NT)
            return;
        end

        if isfield(F,'Het') && ~isempty(F.Het)
            F.Het = aParseHetEntities(F);
        else
            F.Het = [];
        end

        %
        het_ids ={};
        for i = 1:length(F.Het)
            het_ids{end+1}=aGetHetId(F,i);
        end
        if length(unique(het_ids)) ~= length(het_ids)
            F = zAddNTData([pdb_id '.pdb']);
            F.Het = aParseHetEntities(F);
        end
        %

        N = length(F.NT);
        A = length(F.AA);
        H = length(F.Het);
        S = N + A + H;

        c = zeros(S,3);
        for i = 1:N
            if ~isempty(F.NT(i).Center),
                c(i,:) = F.NT(i).Center;
            end
        end

        for i = 1:A
            if ~isempty(F.AA(i).Center),
                c(N+i,:) = F.AA(i).Center;
            end
        end

        for i = 1:H
            if ~isempty(F.Het(i).Center),
                c(N+A+i,:) = F.Het(i).Center;
            end
        end

        F.Distance = zMutualDistance(c,MAXDISTANCE);

        % precompute and store all ids
        ids = cell(1,S);
        for i = 1:S
            ids{i} = aEntityId(i,N,A,F);
        end

        [x,y] = ind2sub(S,find(F.Distance));
%         X = length(x);

        fid = fopen(FILENAME,'w');
        for i = 1:length(x)
            fprintf(fid,'"%s","%s","%.2f"\n',ids{x(i)},ids{y(i)},F.Distance(x(i),y(i)));
        end
        fclose(fid);

        status = 0;

    catch err
        err_msg = sprintf('Error "%s" in aGetDistances on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);
        status = 1;
    end

end

function [id] = aEntityId(i, N, A, F)

    if i <= N
        id = aGetNTId(F,i);
    elseif i <= N+A
        id = aGetAAId(F,i-N);
    else
        id = aGetHetId(F,i-(N+A));
    end

end


%         D = cell(1,X);
%         for i = 1:X
%             D{i}{1} = ids{x(i)};
%             D{i}{2} = ids{y(i)};
%             D{i}{3} = F.Distance(x(i),y(i));
%         end
%         Dist = struct('id1', cell(1,X), ...
%                       'id2', '', ...
%                       'dist', 0);
%
%         for i = 1:X
%             Dist(i).id1  = ids{x(i)};
%             Dist(i).id2  = ids{y(i)};
%             Dist(i).dist = F.Distance(x(i),y(i));
%         end
%
%         L = length(Dist);
%         D = {D};%D(1:20);
