"""This is a module that contains all stages related to loading unit data.

Stages
------
units.info : Populates the `units.info` table.
units.coordinates : Stores the coordinates.
units.distances : Store the distances between units.
units.quality : Fetch and store quality (RSR, RSRZ) data.
units.redundant : Stores information about which units are redundant. Probably
                unneeded
units.centers : Stores the centers or structures.
units.rotation : Stores rotation matrices.
units.loader : Runs all stages here.

"""

