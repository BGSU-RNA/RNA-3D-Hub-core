# PDB files that should be skipped.  Add others as necessary, and note reason
# for exclusion when known.
#
# This is a very large virus file that should be skipped.
# (JJC, 2017-02-09:  no idea what this old comment refers to.  Neither 4V3P
#   nor 4V4G.)
#
# Shorthand notations:
# IPL0 = fails interactions.pairwise; error "Line 0 did not include both units"
#


#       '7PJS' # failing exp_sep.mapping as on 2024-06-20 because it does not have unit_info for chains v,w,z


SKIP = {
      '4V3P'
    , '4V4G'
    , '4ZPY' # failing units.info as of 2017-05-26 (JJC)
    , '5O58' # IPLO 2017-10-25 trial dry run (JJC)
    , '5W3V' # IPLO 2017-10-25 trial dry run (JJC)
    , '6AWB' # IPLO 2017-10-25 trial dry run (JJC)
    , '6AWC' # IPLO 2017-10-25 trial dry run (JJC)
    , '6B0B' # IPLO 2017-10-25 trial dry run (JJC)
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
    , '1DHH' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1DRN' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1FC8' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1GTC' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1HO6' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1HOQ' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '1OKA' # missing pdbx_struct_oper_list block 2020-11-08 (CLZ)
    , '4XK0' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)
    , '5FCI' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)
    , '5JEA' # failing quality.units (missing clashes in validation) as of 2018-01-09 (JJC)
    , '1UTV' # failing loops.extractor as of 2018-01-10 (JJC) # still failing 2019-01-15 (JJC)
    , '5O61' # failing loops.extractor as of 2018-01-10 (JJC) # still failing 2019-01-15 (JJC)
    , '3U5F' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    , '4WUS' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    , '4WWE' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    , '4WWT' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    , '4Z3Q' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)
    , '4Z3R' # failing loops.positions as of 2018-01-10 (JJC) # still failing 2019-01-16 (JJC)


    , '1VQL' # failing quality.units as of 2018-02-09 (JJC)
    , '1VQM' # failing quality.units as of 2018-02-09 (JJC)
    , '1YJW' # failing quality.units as of 2018-02-09 (JJC)
    , '2UU9' # failing quality.units as of 2018-02-09 (JJC)
    , '2UUA' # failing quality.units as of 2018-02-09 (JJC)
    , '2UUB' # failing quality.units as of 2018-02-09 (JJC)
    , '2UUC' # failing quality.units as of 2018-02-09 (JJC)
    , '2VQE' # failing quality.units as of 2018-02-09 (JJC)
    , '2VQF' # failing quality.units as of 2018-02-09 (JJC)
    , '2X2Q' # failing quality.units as of 2018-02-09 (JJC)
    , '2XO0' # failing quality.units as of 2018-02-09 (JJC)
    , '2XZO' # failing quality.units as of 2018-02-09 (JJC)
    , '4Q9Q' # failing quality.units as of 2018-02-09 (JJC)
    , '4Q9R' # failing quality.units as of 2018-02-09 (JJC)
    , '4WZD' # failing quality.units as of 2018-02-09 (JJC)
    , '4YHH' # failing quality.units as of 2018-02-09 (JJC)
    , '5DDQ' # failing quality.units as of 2018-02-09 (JJC)
    , '5FCJ' # failing quality.units as of 2018-02-09 (JJC)
    , '5I8Q' # failing quality.units as of 2018-02-09 (JJC)
    , '5IP2' # failing quality.units as of 2018-02-09 (JJC)
    , '5VP2' # failing quality.units as of 2018-02-09 (JJC)
    , '6C6K' # failing species.mapping as of 2018-05-25 (JJC)
    , '6SAE' # failing export.cifatom as of 2019-12-04 (CLZ)
    , '6WLN' # failing mat_files as of 2020-07-08 (CLZ) Out of memory
    , '6WLO' # failing mat_files as of 2020-07-08 (CLZ) Out of memory

}

#    , '2H0S' # failing loops.quality as of 2018-01-11 (JJC) # works 3/7/2020 CLZ

    # This block of files may work now that exp_seq.mapping is revised to fill in missing chains
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

