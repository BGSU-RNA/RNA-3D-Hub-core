function [MM] = aSymmetrizeMatrix(MM, loop_ids, saveMatFile)

    disp('Making the matrix symmetric...');

    if nargin < 3
        saveMatFile = 1;
    end
    
    N = length(MM);
    Location = getSearchFolder;
    verbose = 1;
        
    BASEPAIR_MISMATCH  = 7;
    BASESTACK_MISMATCH = 8;

    PAIRS  = 1:12;
    STACKS = 21:23;
    SYMMETRIC_PAIRS = [1 2 7 8 11 12]; % cWW, tWW, cHH, tHH, cSS, tSS
        
    FILENAME = 'MM_symmetrize.mat';
    LOGFILE  = 'MM_symmetrize.txt';
    
    fid = fopen(LOGFILE, 'a');

    for iLoop = 1:N
        
        fprintf('%i out of %i\n', iLoop, N);

        ind = find(MM(iLoop,:) <= 2 & MM(iLoop,:) > 0); % less than any disqualified match
        
        for jLoop = ind

            motif1 = loop_ids{iLoop};
            motif2 = loop_ids{jLoop};
            
            switch aCompareInteractions
                case BASEPAIR_MISMATCH
                    MM(iLoop,jLoop) = BASEPAIR_MISMATCH;
                    MM(jLoop,iLoop) = BASEPAIR_MISMATCH; 
                case BASESTACK_MISMATCH
                    MM(iLoop,jLoop) = BASESTACK_MISMATCH;
                    MM(jLoop,iLoop) = BASESTACK_MISMATCH;                
                otherwise                    
                    MM(jLoop,iLoop) = MM(iLoop,jLoop);
            end
            
        end
        
    end
    
    if saveMatFile
        save(FILENAME,'MM','loop_ids');
    end
    fclose(fid);
    checkMatchingMatrix(MM);
    

    function [disqualify] = aCompareInteractions

        disqualify = 0;
        
        filename = getSearchAddress(motif1, motif2);
        if ~exist(filename,'file')        
            filename = getSearchAddress(motif2, motif1);
        end
        load(filename, 'Search');

        cand = find(Search.Discrepancy == min(Search.Discrepancy),1);
        pdb  = Search.Candidates(cand,end);
        nts  = Search.Candidates(cand,1:end-1);
        
        foundEdgesFix = fix(Search.File(pdb).Edge(nts,nts));
        queryEdgesFix = fix(Search.Query.Edge);
        foundEdgesFixAbs = abs(foundEdgesFix);
        queryEdgesFixAbs = abs(queryEdgesFix);

        E = length(foundEdgesFixAbs);

        for i = 1:E

            for j = (i+1):E

                isPair = foundEdgesFixAbs(i,j) <= 12 && queryEdgesFixAbs(i,j) <= 12 && ...
                         queryEdgesFixAbs(i,j) > 0   && foundEdgesFixAbs(i,j) > 0;
                
                if isPair
                    
                    % cWw vs cwW OK, tHS vs tSH not OK
                    is_pair_incompatible = queryEdgesFixAbs(i,j) ~= foundEdgesFixAbs(i,j);
                    
                    % cases like -9 and 9 vs 9 and 9                    
                    is_pair_asymmetric_incompatible = isempty(intersect(foundEdgesFixAbs(i,j), SYMMETRIC_PAIRS)) && ...
                                                      foundEdgesFix(i,j) ~= queryEdgesFix(i,j) && ...
                                                      queryEdgesFixAbs(i,j) == foundEdgesFixAbs(i,j);
                    
                    if is_pair_incompatible || is_pair_asymmetric_incompatible
                        disqualify = BASEPAIR_MISMATCH;
                        if verbose
                            fprintf('%s %s\n',zEdgeText(foundEdgesFixAbs(i,j)), ...
                                              zEdgeText(queryEdgesFixAbs(i,j)));
                        end
                        
                        annotate_conflicting_interactions();
                        
                        return;                        
                    end
                end

                % a coplanar pair in one structure matches a stack in another
                is_pair_stack_conflict = (ismember(foundEdgesFixAbs(i,j), PAIRS)   && ...
                                          ismember(queryEdgesFixAbs(i,j), STACKS)) || ...
                                         (ismember(queryEdgesFixAbs(i,j), PAIRS)   && ...
                                          ismember(foundEdgesFixAbs(i,j), STACKS));

                % disqualify only if the basepair IS coplanar
                if is_pair_stack_conflict
                    % one of the two values is a pair
                    if ismember(foundEdgesFixAbs(i,j), PAIRS)
                        % use Search.Candidates ordering to get the right nucleotides
                        pair = zLooseCoplanar(Search.File(pdb).NT(nts(i)), Search.File(pdb).NT(nts(j)));
                        if isfield(pair, 'Coplanar') && pair.Coplanar ~= 0
                            disqualify = BASESTACK_MISMATCH;
                        end
                    else
                        % queryEdges is in the same order as Search.Query.NT
                        pair = zLooseCoplanar(Search.Query.NT(i), Search.Query.NT(j));
                        if isfield(pair, 'Coplanar') && pair.Coplanar ~= 0
                            disqualify = BASESTACK_MISMATCH;
                        end
                    end
                    
                    if disqualify == BASESTACK_MISMATCH
                        if verbose
                            fprintf('%s %s\n',zEdgeText(foundEdgesFixAbs(i,j)), ...
                                              zEdgeText(queryEdgesFixAbs(i,j)));
                        end

                        annotate_conflicting_interactions();
                    
                        return;
                    end
                end

            end
        end                        
        
        function [] = annotate_conflicting_interactions()
            nt1 = [Search.File(pdb).NT(nts(i)).Base Search.File(pdb).NT(nts(i)).Number];
            nt2 = [Search.File(pdb).NT(nts(j)).Base Search.File(pdb).NT(nts(j)).Number];
            int1 = strtrim(zEdgeText(foundEdgesFixAbs(i,j)));
            nt3 = [Search.Query.NT(i).Base Search.Query.NT(i).Number];
            nt4 = [Search.Query.NT(j).Base Search.Query.NT(j).Number];
            int2 = strtrim(zEdgeText(queryEdgesFixAbs(i,j)));
            message = sprintf('%s %s %s vs %s %s %s', nt1, int1, nt2, nt3, int2, nt4);
            fprintf(fid, '"%s","%s","%i","%s"\n', motif1, motif2, disqualify, message);                        
        end
        
    end

end