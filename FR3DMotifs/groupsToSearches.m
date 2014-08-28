% Form consensus structural alignment for motif groups based on pairwise
% searches.

function [] = groupsToSearches(Location,groups)
    
    destination = [Location filesep 'Groups'];
    if ~exist(destination,'dir')
        mkdir(destination);
    end
    notfound = 0;    
    
    save([Location filesep 'groups.mat'], 'groups');

    for i = 1:length(groups)   
                                  
        disp(i);

        clear Group Search;

        Group.Candidates = [];
        Group.File = [];
        Group.Discrepancy = [];
                                
        for j = 1:length(groups{i})
            
            search1 = getSearchAddress(groups{i}{1}, groups{i}{j});
            search2 = getSearchAddress(groups{i}{j}, groups{i}{1});
            load( getPrecomputedDataAddress(groups{i}{j}) );
                        
            if exist(search1,'file')
                
                load(search1);
                
                if j == 1
                    
                    bulges  = aDetectBulgedBases(Search.File);
                    indices = setdiff(Search.Candidates(1,1:end-1), bulges);
                    Group.Candidates = [indices 1];                   
                    
                else

                    [a,b,c] = intersect(Search.Query.Indices, Group.Candidates(1,1:end-1));
                    b = reshape(b, 1, []);
                    c = reshape(c, 1, []);
                    
                    indices = Search.Candidates(1,b);                        
                    Group.Candidates(end+1,c) = indices;
                    Group.Candidates(end,end) = j;

                end
                
                Group.File = [Group.File File];
                Group.Discrepancy(end+1) = Search.Discrepancy(1);
            
            elseif exist(search2,'file')

                load(search2);

                [a,b,c]=intersect(Group.Candidates(1,1:end-1),Search.Candidates(1,1:end-1));
                b = reshape(b, 1, []);
                c = reshape(c, 1, []);
                
                Group.Candidates(end+1,b) = Search.Query.Indices(c);
                                                                                                
                Group.Candidates(end,end) = j;
                Group.File = [Group.File File];
                Group.Discrepancy(end+1) = Search.Discrepancy(1);                                                
                
            else
                fprintf('Warning! Not found %s, %s\n',groups{i}{1},groups{i}{j});
                notfound = 1;
%                 keyboard;                
                break;
            end            

        end    
        
        if notfound ~= 1
            % delete cols with zeros
            G = length(Group.Candidates(1,:));
            to_delete = zeros(1,G);
            for j = 1:G
                if any(Group.Candidates(:,j)==0)
                    to_delete(j) = j;
                end
            end
            to_delete(to_delete==0) = [];
            Group.Candidates(:,to_delete) = [];
            Group.Query.RelCutoff  = max(Group.Discrepancy);
            Group.Query.DiscCutoff = max(Group.Discrepancy);
            Group.Query.Geometric = 1;
            Group.Query.NumNT = length(Group.Candidates(1,:)) - 1;
            Group.Query.Name  = ['Group_' sprintf('%03d',i)];

            % order by similarity
            if length(Group.Candidates(:,1)) > 1
                Group.DiscComputed = zeros(1,length(Group.Candidates(:,1)));

                Group = xMutualDiscrepancy(Group.File,Group,400);
                p = zOrderbySimilarity(Group.Disc);    
                Group.Disc = Group.Disc(p,p);
                Group.Candidates = Group.Candidates(p,:);

                temp = {Group.File.Filename};
                Group.LoopsOrderedBySimilarity = temp(p);

            else
                Group.Disc = 0;
                Group.LoopsOrderedBySimilarity = {Group.File.Filename};
            end

            % order relative to the exemplar
            exemplar = findExemplar(Group.Disc);
            [Group.Discrepancy, exemplar_order] = sort(Group.Disc(exemplar, :));
            Group.LoopsOrderedByDiscrepancy = Group.LoopsOrderedBySimilarity(exemplar_order);   
            Group.Candidates = Group.Candidates(exemplar_order, :);

            Search = Group;
            
            if ~isfield(Search.File(1),'LooseCoplanar')     
                Search = aAddLooseCoplanar(Search);
            end
            
            Search.consensusEdge = pConsensusInteractions(Search);
            Search.Signature = zMotifSignature(Search.consensusEdge);
            fprintf('%s\n',Search.Signature);                                    
            
            save([destination filesep Group.Query.Name '.mat'], 'Search');
        else
            notfound = 0;
        end
        
    end
        
end

function [index] = findExemplar(M)

    N = length(M);

    vals = zeros(1, N);

    for i = 1:N
        vals(i) = sum(M(i,:));
    end

    [minVal, index] = min(vals);

end

function [Search] = aAddLooseCoplanar(Search)

    [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
    N = N - 1;                              % number of nucleotides

    for ff = 1:length(Search.File),
      F = Search.File(ff);
      if ~isempty(F.NT),
        F.LooseCoplanar = sparse(F.NumNT,F.NumNT);
        NewFile(ff) = F;
      end
    end
    Search.File = NewFile;

    for c = 1:length(Search.Candidates(:,1)),
      ff = Search.Candidates(c,N+1);
      i = Search.Candidates(c,1:N);

      for a = 1:length(i),
        for b = (a+1):length(i),
          if Search.File(ff).Edge(i(a),i(b)) ~= 0,
            NT1 = Search.File(ff).NT(i(a));
            NT2 = Search.File(ff).NT(i(b));
            Pair = zLooseCoplanar(NT1,NT2);
            Search.File(ff).LooseCoplanar(i(a),i(b)) = Pair.Coplanar;
            Search.File(ff).LooseCoplanar(i(b),i(a)) = Pair.Coplanar;
          end
        end
      end
    end

end
