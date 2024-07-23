# PDB files that should be skipped.
# Add others as necessary
# Note reason for exclusion when known.
#
#
# Shorthand notations:
# IPL0 = fails interactions.pairwise; error "Line 0 did not include both units"
#

# Things to continue working on
#       '7PJS' # failing exp_sep.mapping as on 2024-06-20 because it does not have unit_info for chains v,w,z
# try running unit_info on this file to see if it can be fixed

SKIP = {
      '4V3P' # killed 2024-06-22 (CLZ), export.cifatom:Writing cifatom file for 4V3P
    , '4V4G' # Matlab out of memory 2024-06-22 (CLZ)
    , '4ZPY' # exp_seq.info:Nothing to process, virus, maybe no RNA 2024-06-22 (CLZ)
    , '5O58' # IPLO 2017-10-25 trial dry run (JJC) Matlab problem 2024-06-22 (CLZ)
    , '5W3V' # IPLO 2017-10-25 trial dry run (JJC) Matlab problem 2024-06-22 (CLZ)
    , '6AWB' # IPLO 2017-10-25 trial dry run (JJC) Matlab problem 2024-06-22 (CLZ)
    , '6AWC' # IPLO 2017-10-25 trial dry run (JJC) Matlab problem 2024-06-22 (CLZ)
    , '6B0B' # IPLO 2017-10-25 trial dry run (JJC) Matlab problem 2024-06-22 (CLZ)
    # after we replace Matlab with python code, we can try these again
    , '6WLN' # failing mat_files as of 2020-07-08 (CLZ) Out of memory
    , '6WLO' # failing mat_files as of 2020-07-08 (CLZ) Out of memory

    , '5O61' # failing loops.extractor as of 2018-01-10 (JJC) # still failing 2019-01-15 (JJC)
    , '2X2Q' # failing quality.units as of 2018-02-09 (JJC)   # still failing 2024-06-20 (CLZ)

    # The ones below all failed on 2024-06-22, generally something deep about pdbx
    , '1CK5' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1CK8' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1CN8' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1CN9' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1DZS' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1E6T' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1E7X' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1EOR' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1FJF' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1GKV' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1GKW' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1H8J' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1HDW' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1HE0' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1HE6' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1JJM' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1JJN' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1NLE' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1NYZ' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1T42' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1VS9' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '1ZFR' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '2I1C' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '2NVS' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '2NYO' # missing pdbx_struct_oper_list block 2017-11-07 (JJC)
    , '104D' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '169D' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '170D' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1FC8' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1GTC' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1HO6' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1HOQ' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1OKA' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)


    ########## The following pdb files are DNA files. Ding got the premission to skipped the following files on 12/17/2021
    ########## The following files will be fixed and removed in the future.
    , '103D', '105D', '106D', '107D', '108D', '132D', '135D', '136D', '139D', '140D', '141D', '142D'
    , '143D', '148D', '149D', '156D', '171D', '175D', '177D', '179D', '186D', '199D', '1A66', '1A6B'
    , '1A6H', '1A83', '1A8N', '1A8W', '1AC7', '1AC9', '1AF1', '1AFF', '1AFZ', '1AG3', '1AGH', '1AGK'
    , '1AGO', '1AGU', '1AGZ', '1AL9', '1AMD', '1AO1', '1AO9', '1AP1', '1AT4', '1AU6', '1AUL', '1AW4'
    , '1AX6', '1AX7', '1AXL', '1AXO', '1AXP', '1AXU', '1AXV', '1B0S', '1B3P', '1B4Y', '1BBX', '1BCB'
    , '1BCE', '1BDZ', '1BE5', '1BHR', '1BJ6', '1BJD', '1BJH', '1BN9', '1BPS', '1BUB', '1BUF', '1BUT'
    , '1BWG', '1BWT', '1BX5', '1C0Y', '1C11', '1C32', '1C34', '1C35', '1C38', '1C7U', '1C95', '1CFL'
    , '1COC', '1CQO', '1CR3', '1CS2', '1CX3', '1CYZ', '1D18', '1D19', '1D20', '1D3X', '1D42', '1D68'
    , '1D69', '1D6D', '1D70', '1DAU', '1DB6', '1DGO', '1DJD', '1DK6', '1DK9', '1DL4', '1DSA', '1DSI'
    , '1DSM', '1DUF', '1DXA', '1ECU', '1EEG', '1EEK', '1EL2', '1ELN', '1EMQ', '1EO4', '1ESS', '1EU2'
    , '1EU6', '1EVM', '1EVN', '1EVO', '1EW1', '1EZN', '1F3S', '1F4S', '1F5E', '1FJ5', '1FJB', '1FKY'
    , '1FKZ', '1FQP', '1FV8', '1FYI', '1FYY', '1FZL', '1FZS', '1FZX', '1G14', '1G1N', '1G22', '1G4D'
    , '1G5D', '1G5E', '1G5L', '1G7Z', '1G80', '1GCC', '1GIP', '1GIZ', '1GJ0', '1GJ2', '1GN7', '1HM1'
    , '1HRY', '1HRZ', '1HT4', '1HT7', '1HVN', '1HVO', '1HWV', '1HX4', '1HZ0', '1HZ2', '1I34', '1I5V'
    , '1I7V', '1IDX', '1IEK', '1IEY', '1IG4', '1II1', '1IV6', '1J46', '1J5K', '1J5N', '1JDG', '1JJP'
    , '1JO1', '1JRV', '1JRW', '1JS5', '1JS7', '1JU0', '1JUU', '1JVC', '1JVE', '1K1H', '1K1R', '1K29'
    , '1K2J', '1K2K', '1K4X', '1K5E', '1K5F', '1K9H', '1K9L', '1KB1', '1KBD', '1KBM', '1KKV', '1KKW'
    , '1KR8', '1KSE', '1KVH', '1KXS', '1L0R', '1LA8', '1LAE', '1LAI', '1LAQ', '1LAS', '1LCC', '1LCD'
    , '1LEJ', '1LVS', '1LWA', '1M6A', '1MSE', '1MSF', '1MXK', '1MYQ', '1N14', '1N17', '1N1K', '1NK2'
    , '1NK3', '1OLD', '1QBY', '1QCH', '1QDF', '1QDH', '1QDI', '1QDK', '1QE7', '1QSK', '1QSX', '1RCS'
    , '1SAA', '1SJK', '1SKP', '1SLS', '1TAN', '1TNE', '1UQA', '1UQB', '1UQC', '1UQD', '1UQE', '1UQF'
    , '1UQG', '1WAN', '1XUE', '1YUI', '1YUJ', '1ZHU', '201D', '202D', '203D', '204D', '214D', '225D'
    , '226D', '229D', '230D', '2ARG', '2DAU', '2EZD', '2EZE', '2EZF', '2EZG', '2GAT', '2HDC', '2KBD'
    , '2NEO', '2STT', '2STW', '3GAT', '3KBD', '3REC', '4GAT', '4KBD', '5GAT'



}

    # These worked on 2024-06-22 (CLZ)  Not sure what changed
    # , '1DHH' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    # , '1DRN' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)

    # These worked on 2024-06-22 (CLZ)  Not sure what changed
    # , '1UTV' # failing loops.extractor as of 2018-01-10 (JJC) # still failing 2019-01-15 (JJC)
    # , '3U5F' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    # , '4WUS' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    # , '4WWE' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    # , '4WWT' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    # , '4Z3Q' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    # , '4Z3R' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)


    # working as of 2024-06-20 (CLZ)
    # , '1VQL' # failing quality.units as of 2018-02-09 (JJC)
    # , '1VQM' # failing quality.units as of 2018-02-09 (JJC)
    # , '1YJW' # failing quality.units as of 2018-02-09 (JJC)
    # , '2UU9' # failing quality.units as of 2018-02-09 (JJC)
    # , '2UUA' # failing quality.units as of 2018-02-09 (JJC)
    # , '2UUB' # failing quality.units as of 2018-02-09 (JJC)
    # , '2UUC' # failing quality.units as of 2018-02-09 (JJC)
    # , '2VQE' # failing quality.units as of 2018-02-09 (JJC)
    # , '2VQF' # failing quality.units as of 2018-02-09 (JJC)
    # , '2XO0' # failing quality.units as of 2018-02-09 (JJC)
    # , '2XZO' # failing quality.units as of 2018-02-09 (JJC)
    # , '4Q9Q' # failing quality.units as of 2018-02-09 (JJC)
    # , '4Q9R' # failing quality.units as of 2018-02-09 (JJC)
    # , '4WZD' # failing quality.units as of 2018-02-09 (JJC)
    # , '4YHH' # failing quality.units as of 2018-02-09 (JJC)
    # , '5DDQ' # failing quality.units as of 2018-02-09 (JJC)
    # , '5FCJ' # failing quality.units as of 2018-02-09 (JJC)
    # , '5I8Q' # failing quality.units as of 2018-02-09 (JJC)
    # , '5IP2' # failing quality.units as of 2018-02-09 (JJC)
    # , '5VP2' # failing quality.units as of 2018-02-09 (JJC)
    # , '4XK0' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)
    # , '5FCI' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)
    # , '5JEA' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)


    # , '6SAE' # failing export.cifatom  as of 2019-12-04 (CLZ) works 2024-06-20 CLZ
    # , '6C6K' # failing species.mapping as of 2018-05-25 (JJC) works 2024-06-20 CLZ

    # , '2H0S' # failing loops.quality as of 2018-01-11 (JJC) # works 3/7/2020 CLZ

    # This block of files work now 2024-06-20 that exp_seq.mapping is revised to fill in missing chains
    # , '2H0Z' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3BO4' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain D (CLZ)
    # , '3CQS' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3CR1' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3I2Q' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain A (CLZ)
    # , '3I2R' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain A (CLZ)
    # , '3I2S' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain A (CLZ)
    # , '3I2U' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain A (CLZ)
    # , '3IIN' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain C (CLZ)
    # , '3OLA' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3ZD4' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3ZD5' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3ZP8' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '3ZVO' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain V (CLZ)
    # , '4K4Y' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain B (CLZ)
    # , '5IT7' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain 2 (CLZ)
    # , '5J4D' # failing loops.quality as of 2018-01-11 (JJC) # no exp seq pos for chain A (CLZ)

    #    , '1QZA' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '1QZB' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '1T1O' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '2AGN' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '3DG0' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '3DG2' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '3DG4' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '3DG5' # failing units.distances as of 2017-03-18 (JJC) # still failing 2018-01-23 (JJC)
    #    , '3J64' # failing units.distances as of 2017-11-27 (JJC) # still failing 2018-01-23 (JJC)
    #    , '6EM3' # failing units.distances (No data saved) as of 2018-02-28 (JJC)
    #    , '6EM4' # failing units.distances (No data saved) as of 2018-02-28 (JJC)
    #    , '6EM5' # failing units.distances (No data saved) as of 2018-02-28 (JJC)
    #    , '5ZZM' # failing units.distances (No data saved) as of 2018-07-02 (JJC)

    #, '5Z56' # failing ife.info as of 2018-10-03 (JJC)
    #, '5Z57' # failing ife.info as of 2018-10-03 (JJC)
    #, '5Z58' # failing ife.info as of 2018-10-03 (JJC)
    #, '6DZK' # failing ife.info as of 2018-10-03 (JJC)
