% output:
% result.loop_ids is a cell array of annotated loops
% result.best_matches is a cell array of best matching motif exemplars.
% Contains -1 for loops without matches.


function [result] = annotatePdbWithMotifs(pdb_id, motif_type)
%function [] = annotatePdbWithMotifs(pdb_id, exemplar_ids, motif_type)

    exemplar_ids = {'IL_3U5F_089','IL_1S72_070','IL_3NKB_002','IL_1X8W_005','IL_1F7Y_003','IL_1DUH_005','IL_1X8W_004','IL_2ZJR_012','IL_2XZM_059','IL_4A1B_132','IL_3UXR_049','IL_1FJG_022','IL_3UXR_092','IL_3U5H_064','IL_1JBR_002','IL_1S72_068','IL_3BNQ_008','IL_2AW7_030','IL_3UXR_008','IL_3U5H_066','IL_2QBG_072','IL_1FJG_041','IL_2ZJR_081','IL_1S72_072','IL_2QBG_098','IL_2XZM_028','IL_1NBS_008','IL_1NTB_003','IL_2PXV_002','IL_2OZB_002','IL_1NUV_001','IL_2AW7_044','IL_3Q1Q_006','IL_2ZJR_026','IL_1DUQ_005','IL_2ZI0_001','IL_3UXR_111','IL_2ZJR_038','IL_2XZM_053','IL_1MJI_004','IL_3U5H_035','IL_3U5H_029','IL_4A1B_054','IL_3D2V_006','IL_3UXR_027','IL_2OZB_001','IL_3U5H_016','IL_2QBG_093','IL_2VPL_004','IL_2AW7_061','IL_3U5H_086','IL_3U5H_015','IL_3UXR_106','IL_3UXR_061','IL_1S72_053','IL_1S72_069','IL_413D_001','IL_3P59_001','IL_1ET4_009','IL_2GDI_005','IL_4A1B_004','IL_2ZJR_101','IL_3Q3Z_004','IL_2XZM_083','IL_2XZM_031','IL_3U5H_113','IL_3UXR_102','IL_3U5F_017','IL_3UXR_067','IL_3UXR_054','IL_3UXR_039','IL_2QBG_061','IL_3UXR_071','IL_2NZ4_013','IL_2QBG_057','IL_1S72_057','IL_1FJG_063','IL_2XZM_006','IL_2QWY_005','IL_1FJG_026','IL_2Y9B_006','IL_1FJG_013','IL_1FJG_050','IL_3RW6_003','IL_2QBG_077','IL_1S72_037','IL_2HOJ_002','IL_3UXR_025','IL_2AW7_041','IL_2XZM_010','IL_2AW7_037','IL_2AW7_003','IL_3U5H_088','IL_3U5H_081','IL_3U5H_059','IL_3UXR_108','IL_3UXR_101','IL_1M5O_005','IL_2AW7_038','IL_1J2B_001','IL_1FJG_018','IL_1U6B_009','IL_1GID_011','IL_4A1B_029','IL_1S72_017','IL_1S72_044','IL_2ZJR_017','IL_3UXR_076','IL_3UXR_082','IL_4A1B_057','IL_4A1B_070','IL_4A1B_074','IL_2XZM_046','IL_283D_001','IL_3DIR_003','IL_1FJG_017','IL_1FJG_049','IL_1U6B_010','IL_1NBS_003','IL_1LNG_003','IL_3GX5_003','IL_1YKV_001','IL_3BNL_003','IL_1UN6_002','IL_1OOA_002','IL_3OWZ_004','IL_1DUH_003','IL_3U5F_025','IL_3U5F_011','IL_3U5F_054','IL_3U5F_075','IL_1S72_100','IL_2ZJR_058','IL_2ZJR_103','IL_2QBG_009','IL_2QBG_015','IL_2QBG_053','IL_4A1B_058','IL_4A1B_078','IL_4A1B_085','IL_4A1B_099','IL_3Q3Z_001','IL_1VBY_001','IL_3D2V_002','IL_3U5F_007','IL_3Q1Q_002','IL_3Q1Q_007','IL_1FJG_009','IL_1FJG_011','IL_1FJG_019','IL_1FJG_031','IL_1FJG_032','IL_1FJG_037','IL_1FJG_045','IL_1FJG_048','IL_3G78_001','IL_3G78_003','IL_3G78_006','IL_3G78_009','IL_3G78_011','IL_1U6B_008','IL_1U6B_011','IL_1X8W_008','IL_1Y0Q_001','IL_1Y0Q_002','IL_2OIU_004','IL_2Y9B_002','IL_2Y9B_004','IL_255D_001','IL_1MZP_003','IL_4A1C_004','IL_2QBZ_001','IL_2QBZ_002','IL_2QBZ_003','IL_1GID_004','IL_1MFQ_006','IL_3NPB_005','IL_3NPB_006','IL_3KTW_003','IL_3SD3_001','IL_3P59_003','IL_3U5F_021','IL_3U5F_037','IL_3U5F_039','IL_3U5F_042','IL_3U5F_064','IL_3U5F_073','IL_3U5F_077','IL_1S72_011','IL_1S72_013','IL_1S72_015','IL_1S72_036','IL_1S72_054','IL_1S72_073','IL_1S72_093','IL_1S72_097','IL_2ZJR_008','IL_2ZJR_037','IL_2ZJR_073','IL_2QBG_014','IL_2QBG_028','IL_2QBG_056','IL_2QBG_067','IL_2QBG_079','IL_2QBG_100','IL_2QBG_101','IL_2QBG_106','IL_2QBG_108','IL_3UXR_040','IL_3UXR_060','IL_3UXR_074','IL_4A1B_002','IL_4A1B_016','IL_4A1B_033','IL_4A1B_043','IL_4A1B_045','IL_4A1B_060','IL_4A1B_068','IL_4A1B_097','IL_4A1B_119','IL_3U5H_005','IL_3U5H_006','IL_3U5H_011','IL_3U5H_020','IL_3U5H_021','IL_3U5H_023','IL_3U5H_038','IL_3U5H_040','IL_3U5H_053','IL_3U5H_056','IL_3U5H_060','IL_3U5H_068','IL_3U5H_074','IL_3U5H_079','IL_3U5H_080','IL_3U5H_095','IL_3U5H_130','IL_3U5H_132','IL_3U5H_147','IL_3U5H_150','IL_2AW7_009','IL_2AW7_015','IL_2AW7_021','IL_2AW7_022','IL_2AW7_023','IL_2AW7_033','IL_2AW7_053','IL_2AW7_055','IL_2AW7_056','IL_2XZM_019','IL_2XZM_022','IL_2XZM_026','IL_2XZM_047','IL_2XZM_063','IL_2XZM_067','IL_2XZM_077','IL_2XZM_087','IL_3R4F_001','IL_3SIV_002','IL_1NTB_002'};
    exemplar_ids = exemplar_ids(1:10);

    % get all loops from pdb_id
    location = fullfile(getenv('MA_root'), 'PrecomputedData', pdb_id, [motif_type '*.mat']);
    files = dir(location);

    N = length(files);
    E = length(exemplar_ids);
    loop_ids = cell(1, N);

    for i = 1:N
        loop_ids{i} = files(i).name(1:end-4);
    end


    % todo: avoid reloading MM_exemplars for each PDB
    MM_exemplars = aCreateMM(exemplar_ids);
    MM_exemplars = aAnalyzeExtraNucleotides(MM_exemplars, exemplar_ids);
    MM_exemplars = aSymmetrizeMatrix(MM_exemplars, exemplar_ids);

    result = struct;
    result.loop_ids = cell(1, N);
    result.best_matches = cell(1, N);

    for i = 1:N

        % create new matching matrix
        MM = MM_exemplars;

        % fill in the row
        for j = 1:E;
            MM(E+1, j) = pairwiseSearch(loop_ids{i}, exemplar_ids{j});
        end

        % fill in the column
        for j = 1:E;
            MM(j, E+1) = pairwiseSearch(exemplar_ids{j}, loop_ids{i});
        end

        MM(end, end) = 0;

        compared_ids = [exemplar_ids loop_ids{i}];
        MM = aAnalyzeExtraNucleotides(MM, compared_ids);
        MM = aSymmetrizeMatrix(MM, compared_ids);

        % NB! text files with annotations get overwritten every time

        matches = find(MM(end, :)>0 & MM(end,:)<=1);
        if isempty(matches)
            result.best_matches{i} = -1;
        else
            [best_match, index] = min(MM(end, matches));
            result.best_matches{i} = exemplar_ids{matches(index)};
        end

        result.loop_ids{i} = loop_ids{i};

%         keyboard;

    end

    %todo: check no_candidates.txt in pairwiseSearch.m

    keyboard;





end