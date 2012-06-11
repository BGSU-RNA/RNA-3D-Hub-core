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
    
    for i = 1:length(list)

        clear Search
        load( fullfile(groups_location, list{i}) );
        if ~exist('Search','var')
            continue
        end
        fprintf([list{i}, '\n']);
        
        motif_id = list{i}(1:end-4);        
        LoopIds = {Search.File.Filename};
        
        N = length(Search.Candidates(:,1));        
        if ~isfield(Search,'LoopsOrdered') && N > 1
            Search.DiscComputed = zeros(1,N);
            Search = xMutualDiscrepancy(Search.File,Search,300);
            p = zOrderbySimilarity(Search.Disc);    
            Search.Disc = Search.Disc(p,p);
            Search.LoopsOrdered = LoopIds(p);
            save( fullfile(groups_location, list{i}),'Search');  
        end
        if N == 1            
            Search.LoopsOrdered = LoopIds;            
        end
                
        for j = 1:N

            loopid = LoopIds{Search.Candidates(j,end)};

            fprintf(MotifList,'"%s","%s"\n', loopid, motif_id);
                        
            for pos = 1:length(Search.Candidates(j,1:end-1))
                NTid = aGetNTId(Search.File(Search.Candidates(j,end)),Search.Candidates(j,pos));
                fprintf(MotifPositions,'%s,%s,%s,%i\n',motif_id,loopid,NTid,pos);
            end
        
            reordered = find(ismember(Search.LoopsOrdered,loopid)==1,1);

            fprintf(CandOrder,'%s,%s,%i,%i\n',motif_id,loopid,j,reordered);
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