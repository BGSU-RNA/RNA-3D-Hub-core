% Returns a cell array with loop ids making cliques. The first entry in 
% each clique is the motif exemplar, which has the lowest discrepacy
% to all other members of the clique.

function [groups] = aMaximumCliques(M, names, cutoff)

    M_original = M;

    M( M <= cutoff ) = 0;
    M( M >  cutoff ) = 10;

    N = length(M);
    loops = 1:N;
    
    groups = cell(1,N);
    counter = 1;
    
    tempGroups = cell(1,2);

    while ~isempty( find(M==0,1) )

        % clique = rGetMaximalClique(M); % greedy, fast, less accurate
        clique = rCliqueBB(M); % branch and bound, slower, more accurate
        
        if ~all(M(clique,clique)<=cutoff);
            disp('Problem');
        else
            fprintf('Valid clique,%i\n',length(clique));
        end
                
        aCompareAlternativeCliques;
        
        % put exemplar first
        exemplar = findExemplar(M_original(clique, clique));
        temp = clique(1);
        clique(1) = clique(exemplar);
        clique(exemplar) = temp;
        
        M(clique,:) = 10;
        M(:,clique) = 10;
        
        groups{counter} = names(clique);
    %     groups{counter} = clique;
        counter = counter + 1;
    end

    groups(cellfun(@isempty,groups)) = [];    
    
    
    function aCompareAlternativeCliques

            C = length(clique);
            if C == 1
                return;
            end
            MAX_SIZE = C;
            MIN_DISC = sum(sum(M_original(clique,clique))) / C;
            alternative_cliques = {};

            for i = 1:C

                tempM = M;
                tempM(clique(i),clique) = 10;
                tempM(clique,clique(i)) = 10;
                tempM(clique(i),clique(i)) = 0;

                new_clique = rCliqueBB(tempM);
                if length(new_clique) == MAX_SIZE
                    alternative_cliques{end+1} = new_clique;                
                end

            end

            fprintf('\tFound %i alternative cliques\n',length(alternative_cliques));

            for i = 1:length(alternative_cliques)

                new_disc = sum(sum(M_original(alternative_cliques{i},alternative_cliques{i}))) / C;

                if new_disc < MIN_DISC

                    clique = alternative_cliques{i};

                    fprintf('\tSwitching to an alternative clique (old disc %f vs new disc %f)\n',MIN_DISC,new_disc);

                    MIN_DISC = new_disc;                

                end

            end

    end

end

function [index] = findExemplar(M)
    
    N = length(M);
    
    vals = zeros(1, N);
    
    for i = 1:N
        vals(i) = sum(M(i,:)) + sum(M(:,i));
    end

    [minVal, index] = min(vals);

end
