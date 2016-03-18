"""String to use when joining to IFE ids to build a new ID"""
IFE_SEPERATOR = '+'

"""Min number of base pairs for a chain to be structured."""
STRUCTURED_BP_COUNT = 5

"""Min number of bp/nt for a chain to be structured."""
STRUCTURED_BP_PER_NT = 0.5

"""Min ration of external/internal cWW to join structured IFEs"""
IFE_EXTERNAL_INTERNAL_FRACTION = 0.6

"""PDBS used when bootstrapping the database"""
BOOTSTRAPPING_PDBS = (
    "124D",
    "157D",
    "1A34",
    "1DUH",
    "1E4P",
    "1EIY",
    "1EKD",
    "1ET4",
    "1F5H",
    "1F5U",
    "1FCW",
    "1FEU",
    "1FG0",
    "1FJG",
    "1G59",
    "1GID",
    "1GRZ",
    "1J5E",
    "1KOG",
    "1MDG",
    "1S72",
    "1VY4",
    "1X8W",
    "2HOJ",
    "2HOK",
    "2HOL",
    "2HOM",
    "2HOO",
    "2QQP",
    "3CPW",
    "3GPQ",
    "3T4B",
    "4MGM",
    "4MGN",
    "4OAU",
    "4PMI",
    "4V42",
    "4V4Q",
    "4V6M",
    "4V6R",
    "4V7R",
    "4V7W",
    "4V88",
    "4V8G",
    "4V8I",
    "4V9K",
    "4V9O",
    "4V9Q",
)

"""Max allowed discrepancy between chains for NR grouping"""
NR_DISCREPANCY_CUTOFF = 0.5

"""Min percent increase of length to change representative"""
NR_LENGTH_PERCENT_INCREASE = 1.0

"""Min percent increase of bp count to change representative"""
NR_BP_PERCENT_INCREASE = 1.0

"""Length cutoff before being matched as small"""
CORRESPONDENCE_SMALL_CUTOFF = 36

"""Length cutoff before being matched as huge"""
CORRESPONDENCE_HUGE_CUTOFF = 2000
