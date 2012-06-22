function [MM] = aSymmetrizeMatrix(MM, loop_ids)

    disp('Making the matrix symmetric...');

    N = length(MM);
    Location = getSearchFolder;
    verbose = 1;
        
    BASEPAIR_MISMATCH  = 7;
    BASESTACK_MISMATCH = 8;
    PAIRNPAIR_MISMATCH = 9;    

    PAIRS  = 1:12;
    STACKS = 21:23;
    SYMMETRIC_PAIRS = [1 2 7 8 11 12]; % cWW, tWW, cHH, tHH, cSS, tSS
        
    FILENAME = 'MM_symmetrize.mat';
    LOGFILE  = 'MM_symmetrize.txt';
    
    fid = fopen(LOGFILE, 'w');

    for iLoop = 1:N
        
        fprintf('%i out of %i\n', iLoop, N);

        ind = find(MM(iLoop,:) <= 2 & MM(iLoop,:) > 0); % less than any disqualified match
        
        for jLoop = ind

            motif1 = loop_ids{iLoop};
            motif2 = loop_ids{jLoop};
            
%             if isequal(iLoop, jLoop) %|| isequal(MM(iLoop,jLoop), MM(jLoop,iLoop)) %|| MM(iLoop,jLoop) < 0
%                 continue;
%             end                                                   
            
            switch aCompareInteractions
                case BASEPAIR_MISMATCH
                    MM(iLoop,jLoop) = BASEPAIR_MISMATCH;
                    MM(jLoop,iLoop) = BASEPAIR_MISMATCH; 
                case BASESTACK_MISMATCH
                    MM(iLoop,jLoop) = BASESTACK_MISMATCH;
                    MM(jLoop,iLoop) = BASESTACK_MISMATCH;                
%                 case PAIRNPAIR_MISMATCH
%                 MM(iLoop,jLoop) = PAIRNPAIR_MISMATCH;
%                 MM(jLoop,iLoop) = PAIRNPAIR_MISMATCH;                
                otherwise                    
                    MM(jLoop,iLoop) = MM(iLoop,jLoop);
            end
            
        end
        
    end
    
    save(FILENAME,'MM','loop_ids');
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
        
%         if ~isfield(Search.Query,'Edge')
%             load(getPrecomputedDataAddress(Search.Query.Filename));
%             Search.Query.Edge = File.Edge(Search.Query.Indices, Search.Query.Indices);
%             save(filename,'Search');
%         end
        
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

                % pair in one structure matches a stack in another
                is_pair_stack_conflict = ismember(foundEdgesFixAbs(i,j), PAIRS) && ...
                                         ismember(queryEdgesFixAbs(i,j), STACKS);
                                     
                if is_pair_stack_conflict
                    disqualify = BASESTACK_MISMATCH;
                    if verbose                    
                        fprintf('%s %s\n',zEdgeText(foundEdgesFixAbs(i,j)), ...
                                          zEdgeText(queryEdgesFixAbs(i,j)));
                    end
                    
                    annotate_conflicting_interactions();                    
                    
                    return;
                end


    %             if (queryEdgesFixAbs(i,j) <= 12 && queryEdgesFixAbs(i,j) > 0) || (foundEdgesFixAbs(i,j) <= 12 && foundEdgesFixAbs(i,j) > 0)
    %                 
    %                 if foundEdgesFixAbs(i,j) ~= queryEdgesFixAbs(i,j) 
    %                     if foundEdgesFixAbs(i,j) ~= queryEdgesFixAbs(i,j)+100 
    %                         if foundEdgesFixAbs(i,j)+100 ~= queryEdgesFixAbs(i,j)
    % %                             fprintf('%s %s\n',zEdgeText(foundEdgesFixAbs(i,j)),zEdgeText(queryEdgesFixAbs(i,j)));    
    %                             disqualify = PAIRNPAIR_MISMATCH;
    %                             return;
    %                         end
    %                     end
    %                 end
    %                 
    %             end



            end
        end                        
        
        function [] = annotate_conflicting_interactions()
            nt1 = [Search.File(pdb).NT(i).Base Search.File(pdb).NT(i).Number];
            nt2 = [Search.File(pdb).NT(j).Base Search.File(pdb).NT(j).Number];
            int1 = strtrim(zEdgeText(foundEdgesFixAbs(i,j)));
            nt3 = [Search.Query.NT(i).Base Search.Query.NT(i).Number];
            nt4 = [Search.Query.NT(j).Base Search.Query.NT(j).Number];
            int2 = strtrim(zEdgeText(queryEdgesFixAbs(i,j)));
            message = sprintf('%s %s %s vs %s %s %s', nt1, int1, nt2, nt3, int2, nt4);
            fprintf(fid, '"%s","%s","%i","%s"\n', motif1, motif2, disqualify, message);                        
        end
        
    end

end