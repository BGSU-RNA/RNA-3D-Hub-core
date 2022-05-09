function [] = exportMotifRelease(location)

    groups_location = fullfile(location, 'Groups');

    temp = dir( fullfile(groups_location, '*.mat') );
    list = {temp.name};
    if isempty(list)
        fprintf('Check your input\n');
        return;
    end

    MotifList      = fopen([location filesep 'MotifList.csv'],'w');
    MotifPositions = fopen([location filesep 'MotifPositions.csv'],'w');
    CandOrder      = fopen([location filesep 'MotifLoopOrder.csv'],'w');
    DiscValues     = fopen([location filesep 'MutualDiscrepancy.csv'],'w');
    BpSignatures   = fopen([location filesep 'MotifBpSignatures.csv'],'w');

    fprintf('exportMotifRelease: There are %d motif groups\n', length(list))
    num_loops = 0;

    for i = 1:length(list)

        clear Search
        load( fullfile(groups_location, list{i}) );
        if ~exist('Search','var')
            continue
        end
        fprintf([list{i}, '\n']);

        motif_id = list{i}(1:end-4);
        LoopIds = {Search.File.Filename};

        fprintf(BpSignatures, '"%s","%s"\n', motif_id, Search.Signature);

        N = length(Search.Candidates(:,1));

        num_loops = num_loops + N;

        for j = 1:N

            loopid = Search.LoopsOrderedByDiscrepancy{j};

            fprintf(MotifList,'"%s","%s"\n', loopid, motif_id);

            for pos = 1:length(Search.Candidates(j,1:end-1))
                NTid = Search.File(Search.Candidates(j,end)).NT(Search.Candidates(j,pos)).ID;
                fprintf(MotifPositions,'%s,%s,%s,%i\n',motif_id,loopid,NTid,pos);
            end

            reordered = find(ismember(Search.LoopsOrderedBySimilarity,loopid)==1,1);

            fprintf(CandOrder,'%s,%s,%i,%i\n',motif_id,loopid,j,reordered);
            if N > 1
                for k = 1:N
                    fprintf(DiscValues,'%s,%.04f,%s\n',....
                        loopid,Search.Disc(reordered,k),Search.LoopsOrderedBySimilarity{k});
                end
            end

        end

    end

    fprintf('exportMotifRelease: There are %d loops\n', num_loops)

    fclose(MotifList);
    fclose(MotifPositions);
    fclose(CandOrder);
    fclose(DiscValues);
    fclose(BpSignatures);

    fprintf('exportMotifRelease: wrote out all data files for the motif release.\n');

end