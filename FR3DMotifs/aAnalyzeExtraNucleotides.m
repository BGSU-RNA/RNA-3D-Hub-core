function [MM] = aAnalyzeExtraNucleotides(MM, loop_ids, output_dir, saveMatFile)

    disp('aAnalyzeExtraNucleotides: Analyzing extra nucleotides...');

    tic

    if nargin < 4
        saveMatFile = 1;
    end

    Location = getSearchFolder;

    % currently the penalty must be integer because of the database structure
    % penalties are also set in aSymmetrizeMatrix.m
    % these numbers keep them distinct; helps to track down reasons
    HAIRPIN_STACK_PENALTY1 = 3;
    HAIRPIN_STACK_PENALTY2 = 3;
    BP_PENALTY             = 4;
    NEAR_BP_PENALTY        = 5;
    STACK_PENALTY          = 6;
    MISMATCHED_BULGE       = 9;

    FILENAME = fullfile(output_dir,'MM_extraNTs.mat');
    LOGFILE  = fullfile(output_dir,'MM_extraNTs.log');

    % internal codes for basepairs, stacking; only positive values
    BASEPAIRS  = 1:12;
    STACKS     = [21:23 121:123];
    NEAR_PAIRS = 101:112;

    isHairpin = strcmp(loop_ids{1}(1:2), 'HL');
    isInternal = strcmp(loop_ids{1}(1:2), 'IL');

    N = length(MM(1,:));

    fid = fopen(LOGFILE, 'a');

    fprintf(fid,'MM is a %d by %d matrix\n', N, N);

    for i = 1:N

        fprintf(fid,'aAnalyzeExtraNucleotides: Checking loop %s, %i out of %i against all matched loops\n',loop_ids{i},i,N);

        % indices of reasonable matches found
        ind = find( MM(i,:) > 0 & MM(i,:) < HAIRPIN_STACK_PENALTY1 );

        % load loop i
        if length(ind) > 0
            load(getPrecomputedDataAddress(loop_ids{i}), 'File');
            Fi = File; % original coordinates of loop i
            iBulged = aDetectBulgedBases(Fi);
            indicesi = {Fi.NT.Number};
            fprintf(fid,'Loop %s has %d nucleotides and %d bulges\n', loop_ids{i}, length(indicesi), length(iBulged));
        else
            fprintf(fid,'Loop %s has no matching loops\n', loop_ids{i});
        end

        for j = ind

            % load results of searching non-bulged nts of i within j
            load(getSearchAddress(loop_ids{i}, loop_ids{j}), 'Search');

            % find the lowest-discrepancy match
            cand = find(Search.Discrepancy == min(Search.Discrepancy));

            % index of the File with minimum discrepancy
            pdb = Search.Candidates(cand(1),end); %#ok<FNDSB>
            F1  = Search.File(pdb); % nts of loop j found when i is searched in j

            load(getPrecomputedDataAddress(loop_ids{j}), 'File');
            F2 = File; % original coordinates of loop j

            % indices of the nts in j that are matched by nts in i
            coreNts  = Search.Candidates(cand(1),1:end-1);

            % cell array of core nucleotide numbers
            indices1 = {F1.NT(coreNts).Number};

            % cell array of all numbers from the original PDB file
            indices2 = {F2.NT.Number};

            % when there are no extra nts in j compared to i, stop checking
            % recall that bulged nts in i were not part of the search of i within j
            if length(indices1) == length(indices2)
                continue;
            end

            % cell arrays for the chain identifiers for the nts
            chains1 = {F1.NT(coreNts).Chain};
            chains2 = {F2.NT.Chain};

            % cell arrays for model numbers
            if isfield(F1.NT(1),'ModelNum')
                models1 = {F1.NT(coreNts).ModelNum};
                models2 = {F2.NT.ModelNum};
            else
                models1 = {};
                models2 = {};
                models1(1:length(indices1)) = {''};
                models2(1:length(indices2)) = {''};
            end

            indices1 = strcat(models1, indices1, chains1);
            indices2 = strcat(models2, indices2, chains2);

            % extra nucleotides in the original file of loop j
            [extra,indExtra] = setdiff(indices2, indices1);
            indExtra         = reshape(indExtra, 1, []);

            % collect together any interactions made by these extra nucleotides
            F2.Edge = fix(abs(F2.Edge));
            interactions         = reshape(F2.Edge(indExtra,:),1,[]);
            interactionsWithCore = reshape(F2.Edge(indExtra,coreNts),1,[]);
            interactionsInExtra  = reshape(F2.Edge(indExtra,indExtra),1,[]);

            % extra nucleotides make basepairs with any nucleotide in the motif
            if ~isempty(intersect(interactions,BASEPAIRS))

                MM(i,j) = BP_PENALTY;
                MM(j,i) = BP_PENALTY;
                annotate_interaction_conflict(BP_PENALTY, BASEPAIRS);

            % extra nucleotides make near bps with the core nucleotides
            elseif ~isempty(intersect(interactionsWithCore,NEAR_PAIRS))

                MM(i,j) = NEAR_BP_PENALTY;
                MM(j,i) = NEAR_BP_PENALTY;
                annotate_interaction_conflict(NEAR_BP_PENALTY, NEAR_PAIRS);

            % extra nucleotides make more than 1 stacking interaction
            % with the core nucleotides
            elseif length(find(ismember(interactionsWithCore,STACKS))) > 1

                MM(i,j) = STACK_PENALTY;
                MM(j,i) = STACK_PENALTY;
                annotate_interaction_conflict(STACK_PENALTY, STACKS);

            % extra nucleotides in hairpins stack on core nucleotides
            elseif isHairpin && length(find(ismember(interactionsWithCore,STACKS))) >= 1

                MM(i,j) = HAIRPIN_STACK_PENALTY1;
                MM(j,i) = HAIRPIN_STACK_PENALTY1;
                annotate_extra_stacked_nucleotides(HAIRPIN_STACK_PENALTY1);

            % extra nucleotides in a hairpin stack on each other
            elseif isHairpin && length(find(ismember(interactionsInExtra,STACKS))) >= 1

                MM(i,j) = HAIRPIN_STACK_PENALTY2;
                MM(j,i) = HAIRPIN_STACK_PENALTY2;
                annotate_extra_stacked_nucleotides(HAIRPIN_STACK_PENALTY2);

            % if one of these loops has 5 nucleotides
            elseif isInternal && (length(indicesi) == 5 || length(indices2) == 5)

                jBulged = aDetectBulgedBases(F2);

                fprintf(fid,'Considering %s which has %d nts and %d bulged bases\n', loop_ids{j}, length(indices2), length(jBulged));

                if length(indicesi) == 5 && length(iBulged) == 1
                    if length(indices2) == 5 && length(jBulged) == 1 && F2.NT(jBulged).Base ~= Fi.NT(iBulged).Base
                        MM(i,j) = MISMATCHED_BULGE;
                        MM(j,i) = MISMATCHED_BULGE;
                        annotate_mismatched_bulges(MISMATCHED_BULGE);
                    elseif length(indices2) > 5
                        MM(i,j) = MISMATCHED_BULGE;
                        MM(j,i) = MISMATCHED_BULGE;
                        annotate_mismatched_bulges(10);  % code for annotation only
                    elseif length(jBulged) == 0
                        MM(i,j) = MISMATCHED_BULGE;
                        MM(j,i) = MISMATCHED_BULGE;
                        annotate_mismatched_bulges(11);  % code for annotation only
                    end
                elseif length(indices2) == 5 && length(jBulged) == 1
                    if length(indicesi) > 5
                        MM(i,j) = MISMATCHED_BULGE;
                        MM(j,i) = MISMATCHED_BULGE;
                        annotate_mismatched_bulges(10);  % code for annotation only
                    elseif length(iBulged) == 0
                        MM(i,j) = MISMATCHED_BULGE;
                        MM(j,i) = MISMATCHED_BULGE;
                        annotate_mismatched_bulges(11);  % code for annotation only
                    end
                end
            end
        end
    end

    toc

    if saveMatFile
        save(FILENAME, 'MM', 'loop_ids');
    end
    fclose(fid);

    % nested functions

    % output format: NT1 interaction NT2, NT3 interaction NT4
    function [] = annotate_interaction_conflict(penalty, disallowed_interactions)

        comment = {};

        for nt1 = indExtra
            [a,b] = intersect(F2.Edge(nt1,:), disallowed_interactions);
            b = reshape(b, 1, []);
            % extra nucleotides can make more than one basepair
            for nt2 = b
                nt1_name = sprintf('%s%s', F2.NT(nt1).Base, F2.NT(nt1).Number);
                nt2_name = sprintf('%s%s', F2.NT(nt2).Base, F2.NT(nt2).Number);
                bp = strtrim(zEdgeText(F2.Edge(nt1, nt2)));
                comment{end+1} = sprintf('%s %s %s', nt1_name, bp, nt2_name);
            end
        end

        fprintf(fid, '"%s","%s","%i","%s"\n', loop_ids{i}, loop_ids{j}, penalty, aImplode(comment));

    end


    function [] = annotate_extra_stacked_nucleotides(penalty)

        comment = {};

        for nt = indExtra
            comment{end+1} = sprintf('%s%s', F2.NT(nt).Base, F2.NT(nt).Number);
        end

        fprintf(fid, '"%s","%s","%i","%s"\n', loop_ids{i}, loop_ids{j}, penalty, aImplode(comment));

    end


    function [] = annotate_mismatched_bulges(penalty)

        comment = {};

        for iB = iBulged
            comment{end+1} = [Fi.NT(iB).Base Fi.NT(iB).Number];
        end

        comment{end+1} = 'does not match';

        for jB = jBulged
            comment{end+1} = [F2.NT(jB).Base F2.NT(jB).Number];
        end

        fprintf(fid, '"%s","%s","%i","%s"\n', loop_ids{i}, loop_ids{j}, penalty, aImplode(comment));

    end
end
