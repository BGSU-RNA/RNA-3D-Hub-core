"""String to use when joining to IFE ids to build a new ID"""
IFE_SEPERATOR = '+'

"""Min number of base pairs for a chain to be structured."""
STRUCTURED_BP_COUNT = 5

"""Min number of bp/nt for a chain to be structured."""
STRUCTURED_BP_PER_NT = 0.5

"""Min ration of external/internal cWW to join structured IFEs"""
IFE_EXTERNAL_INTERNAL_FRACTION = 0.6

"""PDBS used when bootstrapping the database"""
BOOTSTRAPPING = {
    'args': {
        'skip_stage': ['pdbs.obsolete', 'export.interactions', 'export.loops'],
        'seed': 1,
        'log_level': 'debug',
    },
    'steps': [
        {
            'stage': 'update',
            'args': {},
            'pdbs': ["157D", "1CGM", "1DUH", "1EIY", "1EKD", "1ET4", "1F5H",
                     "1G59", "1GID", "1I9K", "1IBK", "1J5E", "1KOG", "1MDG",
                     "1UTD", "1VY4", "1WMQ", "1X8W", "2G32", "2HOJ", "2HOK",
                     "2HOL", "2HOM", "2HOO", "2IL9", "3CW5", "4A3G", "4A3J",
                     "4CS1", "4FTE", "4PMI", "4Q0B", "4V6R", "4V7R", "4V7W",
                     "4V88", "4V8G", "4V8I", "4V9K"],
        },

        {
            'stage': 'update',
            'args': {},
            'pdbs': ["124D", "157D", "1AQ3", "1CGM", "1DUH", "1E4P", "1EIY",
                     "1EKD", "1ET4", "1F5H", "1G59", "1GID", "1I9K", "1IBK",
                     "1J5E", "1KOG", "1MDG", "1UTD", "1VY4", "1WMQ", "1X8W",
                     "2G32", "2HOJ", "2HOK", "2HOL", "2HOM", "3CW5", "3GPQ",
                     "3J9M", "3T4B", "4A3G", "4A3J", "4CS1", "4FTE", "4PMI",
                     "4Q0B", "4V6R", "4V7R", "4V7W", "4V88", "4V8G", "4V8H",
                     "4V8I", "4V9K", "4V9O", "4V9Q"]
        },

        {
            'stage': 'update',
            'args': {},
            'pdbs': ["157D", "1A34", "1AQ3", "1CGM", "1DUH", "1E4P", "1EIY",
                     "1EKD", "1F5U", "1FCW", "1FEU", "1FG0", "1FJG", "1G59",
                     "1GID", "1GRZ", "1I9K", "1IBK", "1J5E", "1KOG", "1MDG",
                     "1S72", "1UTD", "1VY4", "1WMQ", "1X8W", "2G32", "2HOJ",
                     "2HOK", "2HOL", "2IL9", "2QQP", "2UUA", "3CPW", "3CW5",
                     "3GPQ", "3J9M", "3T4B", "4A3G", "4A3J", "4CS1", "4FTE",
                     "4MCE", "4MGM", "4MGN", "4NGG", "4NMG", "4OAU", "4OQ8",
                     "4OQ9", "4PMI", "4Q0B", "4R3I", "4TUE", "4V42", "4V4Q",
                     "4V6F", "4V6M", "4V6R", "4V7R", "4V7W", "4V88", "4V8G",
                     "4V8H", "4V9O", "4V9Q", "4X4N", "4YBB", "5AJ3"]
        },

        {
            'stage': 'update',
            'args': {},
            'pdbs': ["124D", "157D", "1A34", "1AQ3", "1CGM", "1DUH", "1E4P",
                     "1EIY", "1EKD", "1ET4", "1F5H", "1F5U", "1FCW", "1FEU",
                     "1FG0", "1FJG", "1G59", "1GID", "1GRZ", "1I9K", "1IBK",
                     "1J5E", "1KOG", "1MDG", "1S72", "1UTD", "1VY4", "1WMQ",
                     "1X8W", "1YZ9", "1Z7F", "2G32", "2HOJ", "2HOK", "2HOL",
                     "2HOM", "2HOO", "2IL9", "2MKN", "2QQP", "2UUA", "3CPW",
                     "3CW5", "3GPQ", "3J9M", "3T4B", "4A3G", "4A3J", "4CS1",
                     "4FTE", "4MCE", "4MGM", "4MGN", "4NGG", "4NMG", "4OAU",
                     "4OQ8", "4OQ9", "4PMI", "4Q0B", "4R3I", "4TUE", "4V42",
                     "4V4Q", "4V6F", "4V6M", "4V6R", "4V7R", "4V7W", "4V88",
                     "4V8G", "4V8H", "4V8I", "4V9K", "4V9O", "4V9Q", "4X4N",
                     "4YBB", "5AJ3"]
        }
    ]
}

"""Max allowed discrepancy between chains for NR grouping"""
NR_DISCREPANCY_CUTOFF = 0.4

"""Min percent increase of length to change representative"""
NR_LENGTH_PERCENT_INCREASE = 0.5

"""Min percent increase of bp count to change representative"""
NR_BP_PERCENT_INCREASE = 0.5

"""Length cutoff before being matched as small"""
CORRESPONDENCE_SMALL_CUTOFF = 36

"""Length cutoff before being matched as huge"""
CORRESPONDENCE_HUGE_CUTOFF = 2000

CORRESPONDENCE_EXACT_CUTOFF = 19
CORRESPONDENCE_LIMITED_CHANGES = 80

"""Set of pairs that must be joined by into a equivelance class"""
EQUIVALENT_PAIRS = set([
    (('1S72', '0'), ('1FG0', 'A')),
    (('1S72', '0'), ('1FFZ', 'A')),
])

"""Number of basepairs that an interaction must, ie greater-than, cross to be
long range.
"""
LONG_RANGE = 3

"""NCBI Taxon id for synthetic constructs"""
SYNTHENIC_SPECIES_ID = 32630

"""Min size for groups that must have a homogenous species"""
NR_MIN_HOMOGENEOUS_SIZE = 70

"""Resolution cutoffs to use for NR classes"""
RESOLUTION_GROUPS = [
    '1.5',
    '2.0',
    '2.5',
    '3.0',
    '3.5',
    '4.0',
    '20.0',
    'all'
]

"""The pattern to use for naming nr classes"""
NR_CLASS_NAME = 'NR_{resolution}_{handle}.{version}'

"""The pattern to use for naming motif groups"""
MOTIF_GROUP_NAME = '{type}_{handle}.{version}'

"""Filename to cache to NR data to"""
NR_CACHE_NAME = 'nr'

"""Filename to cache chain_chain comparison data to"""
CCC_CACHE_NAME = 'ccc'

"""Max discrepancy to allow for chain chain discrepancies"""
MAX_RESOLUTION_DISCREPANCY = 20.0

"""Min number of nts to need for computing discrepancies"""
MIN_NT_DISCREPANCY = 3

NR_REPRESENTATIVE_METHOD = 'compscore'
"""What representative selection method to use in nr.builder"""

RSRZ_PAIRED_OUTLIERS = 1
"""The cutoff for finding if a loop contains poorly modeled interactions"""

RSRZ_FICTIONAL_CUTOFF = 1
"""The cutoff for finding if a loop is fictional by RSRZ"""

NR_ALLOWED_METHODS = set(['X-RAY DIFFRACTION'])
"""List of methods are allowed when selecting the representative of an NR
group"""

MOTIF_ALLOWED_METHODS = set(['X-RAY DIFFRACTION'])
#MOTIF_ALLOWED_METHODS = set(['X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY'])
"""list of methods allowed when selected loops from representatives."""

MOTIF_RESOLUTION_CUTOFF = '3.5'
"""The resolution cutoff for representatives to use in the ML atlas"""

MANUAL_IFE_REPRESENTATIVES = set()
#MANUAL_IFE_REPRESENTATIVES = set([
#    '4LFB|1|A',  # TTh SSU
#    '4V9F|1|0',  # Hm LSU
#])
"""
List of IFEs to use as the representative of its equivalence class
Until superseded by some criteria that the pipeline is able to check
"""
# two entries inactivated on 2018-01-09 (JJC)

WORSE_THAN_MANUAL_IFE_REPRESENTATIVES = set()
"""
The collection of automatically selected representatives that have already
been evaulated as worse than the manually selected ones, despite what the
program things.
"""

"""
The coefficients below were set by Craig Zirbel in August 2017 with the
goal of making each quality indicator contribute equally to the standard
deviation of the composite quality score.
"""
COMPSCORE_COEFFICENTS = {
    'resolution': 1,
    'average_rsr': 8,
    'percent_clash': 0.6,
    'average_rscc': 8,
    'rfree': 18,
    'fraction_unobserved': 4,
}

"""
Set this variable to True when you want to rank every single equivalence class and write those rankings to nr_chains, for example, when the ranking system is changed
"""
WRITE_ALL_EQUIVALENCE_CLASS_RANKINGS =  False