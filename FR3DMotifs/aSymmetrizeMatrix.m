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


    for iLoop = 1:N
        
        fprintf('%i out of %i\n', iLoop, N);

        ind = find(MM(iLoop,:) <= 2); % less than any disqualified match
        
        for jLoop = ind

            motif1 = loop_ids{iLoop};
            motif2 = loop_ids{jLoop};
            
            if isequal(iLoop, jLoop) %|| isequal(MM(iLoop,jLoop), MM(jLoop,iLoop)) %|| MM(iLoop,jLoop) < 0
                continue;
            end                                                   
            
            switch aCompareInteractions
                case BASEPAIR_MISMATCH
                    MM(iLoop,jLoop) = BASEPAIR_MISMATCH;
                    MM(jLoop,iLoop) = BASEPAIR_MISMATCH; 
                case BASESTACK_MISMATCH
                    MM(iLoop,jLoop) = BASESTACK_MISMATCH;
                    MM(jLoop,iLoop) = BASESTACK_MISMATCH;                
                case PAIRNPAIR_MISMATCH
%                 MM(iLoop,jLoop) = PAIRNPAIR_MISMATCH;
%                 MM(jLoop,iLoop) = PAIRNPAIR_MISMATCH;                
                otherwise                    
                    MM(jLoop,iLoop) = MM(iLoop,jLoop);
            end
            
        end
        
    end
    
    save(FILENAME,'MM','loop_ids');
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
    end

end