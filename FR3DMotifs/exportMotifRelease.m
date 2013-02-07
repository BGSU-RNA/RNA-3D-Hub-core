function [] = exportMotifRelease(location)
   
    groups_location = fullfile(location, 'Groups');

    temp = dir( fullfile(groups_location, '*.mat') );
    list = {temp.name};
    if isempty(list)
        fprintf('Check your input\n');
        return;
    end
        
    MotifList  = fopen([location filesep 'MotifList.csv'],'w');
    MotifPositions = fopen([location filesep 'MotifPositions.csv'],'w'); 
    CandOrder  = fopen([location filesep 'MotifLoopOrder.csv'],'w'); 
    DiscValues = fopen([location filesep 'MutualDiscrepancy.csv'],'w');
    BpSignatures = fopen([location filesep 'MotifBpSignatures.csv'],'w');
    
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
        if ~isfield(Search,'LoopsOrdered') && N > 1
            Search.DiscComputed = zeros(1,N);
            Search = xMutualDiscrepancy(Search.File,Search,300);
            p = zOrderbySimilarity(Search.Disc);    
            Search.Disc = Search.Disc(p,p);

            % get correct loop_id ordering
            discOrder = reshape(Search.Candidates(:, end), 1, []);
            Search.LoopsOrdered = LoopIds(discOrder);

            % apply similarity ordering to correct loop_ids
            Search.LoopsOrdered = Search.LoopsOrdered(p);
            save( fullfile(groups_location, list{i}),'Search');  
        end
        if N == 1            
            Search.LoopsOrdered = LoopIds;            
            Search.Disc = 0;
        end
                
        % determine final candidate ordering
        % initial ordering is based on discrepancies from the all-against-all searches, 
        % but those values are calculated for all matched nucleotides, not just the core
        % kept in the motif group, so they can be higher than the final values.
        exemplar = find(ismember(Search.LoopsOrdered, Search.File(1).Filename)==1, 1);
        [a, final_order] = sort(Search.Disc(exemplar, :));
        FinalLoopIds = Search.LoopsOrdered(final_order);        
        
        for j = 1:N

            loopid = LoopIds{Search.Candidates(j,end)};

            fprintf(MotifList,'"%s","%s"\n', loopid, motif_id);
                        
            for pos = 1:length(Search.Candidates(j,1:end-1))
                NTid = aGetNTId(Search.File(Search.Candidates(j,end)),Search.Candidates(j,pos));
                fprintf(MotifPositions,'%s,%s,%s,%i\n',motif_id,loopid,NTid,pos);
            end
        
            reordered = find(ismember(Search.LoopsOrdered,loopid)==1,1);
            final_loop_order = find(ismember(FinalLoopIds,loopid)==1,1);

            fprintf(CandOrder,'%s,%s,%i,%i\n',motif_id,loopid,final_loop_order,reordered);
            if N > 1
                for k = 1:N
                    fprintf(DiscValues,'%s,%.04f,%s\n',....
                        loopid,Search.Disc(reordered,k),Search.LoopsOrdered{k});
                end
            end
        
        end
        
    end

    fclose(MotifList);
    fclose(MotifPositions);
    fclose(CandOrder);
    fclose(DiscValues);
    fprintf('Done\n');

end