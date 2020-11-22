"""A module to load unit rotation matrices for modified RNA into the database.
"""

from itertools import groupby
from pymotifs import core
from pymotifs import models as mod
from pymotifs.units.rotation import Loader as RotationLoader


# list the modified nucleotides that should be checked to see if they have a rotation matrix
# if not, one will be computed
# there needs to be a mapping of atoms in each modified nucleotide to the parent atoms,
# in fr3d-python.fr3d.modified_parent_mapping.py
# Apparently this variable needs to be a tuple for the SQL query in get_units_without_rotation
mod_ref = ('DA','DC','DG','DT','A2M', 'CAR', '1MG', '5MC', '5MU', 'H2U', '4SU', 'MIA', 'U34', 'G7M', 'QUO', '1MA', '7MG', 'UDP', 'SRA', 'A23', 'GDP', 'C31', 'U31', 'PGP', 'CCC', 'N', '2MA', 'OMC', '2MG', 'M2G', 'OMG', 'YYG', 'YG', '5BU', 'CG1', 'I', 'T6A', 'PU', 'PPU', '12A', '2MU', '70U', 'U8U', 'PYY', 'IC', 'IG', 'PSU', 'OMU', 'ATL', 'LCA', 'LCC', 'LCG', 'LKC', 'TLN', 'MAD', 'GAO', 'UAR', '5IC', 'FHU', 'CH', 'CBV', 'IU', 'AVC', '5CG', 'P5P', 'GH3', 'GOM', 'SSU', 'PQ1', 'PUY', 'AET', '0C', '0G', '0U', 'UR3', '3DA', 'MTU', 'AP7', '5AA', 'F3N', '31H', '31M', '23G', 'MA6', 'MNU', 'A2L', 'C2L', 'G2L', 'U2L', 'RIA', 'A5M', 'M5M', 'UD5', 'S4C', 'FMU', 'PYO', 'SUR', 'P1P', '6IA', 'RTP', '5FU', 'N6G', 'NF2', 'T2T', 'O2C', 'ONE', 'ZAD', 'ZCY', 'ZGU', 'ZHP', 'ZTH', '6MZ', 'CM0', 'A5O', 'RSP', 'RSQ', 'UPV', 'A44', 'C5L', 'G48', 'T39', '2AU', '10C', '4OC', '8AN', 'ZBC', 'ZBU', 'A3P', 'G46', '1SC', 'N5M', 'US5', '1DP', '3TD', '3AU', 'F3O', '5CF', 'A6A', 'A6C', 'A6G', 'A6U', 'GRB', 'RUS', '574', 'DJF', 'C43', 'U36', 'LMS', 'GDO', 'RPC', 'A9Z', 'UBD', 'U37', '1RN', 'A7E', '7AT', '6FC', '6FU', 'URU', '0A', '0U1', '2SG', 'Y5P', '45A', 'M3O', '50L', '50N', '56B', 'U5M', 'CVC', '8AZ', '6F7', '6NW', '6OO', '6OP', 'F2T', '5HM', '73W', '75B', '8OS', '8RJ', 'G4P', '9QV', 'MUM', 'C4J', 'M1Y', 'A7C', 'EQ4', 'U4M', 'UFB', 'UOA', 'UOB', 'B8H', 'B8K', 'B8Q', 'B8T', 'B8W', 'B9B', 'B9H', 'BGH', 'E6G', 'E7G', 'I4U', 'M7A', 'MHG', 'P4U', 'P7G', '4AC', 'B8N', 'E3C', 'RY', 'KGV', 'JMH', 'KAK', 'LHH', 'LV2', 'UY4', 'UY1', 'O2Z', 'U23')

class Loader(core.SimpleLoader):
    """A class to load rotation matrices into the database.
    """

    # dependencies = set([RotationLoader])  # uncomment this if you want to run on new files, which you shouldn't
    allow_no_data = True                    # pipeline won't crash when we pass back no data
    # has_data = False


    # Organize a list of unit ids into a list of lists, one PDB id at a time
    def group_units(self, units_list):

        return [list(i) for j, i in groupby(units_list, lambda a: a.split('|')[0])]


    def get_units_without_rotation(self,pdb):
        """ Query the database to find all units with no rotation matrix.
        This will only return units from the tuple mod_ref
        """

        with self.session() as session:

            info = mod.UnitInfo
            rotation = mod.UnitRotations

            query = session.query(info.unit_id).\
                            outerjoin(rotation, rotation.unit_id == info.unit_id).\
                            filter(rotation.unit_id == None).\
                            filter(info.unit.in_(mod_ref))

        return [r.unit_id for r in query]


    def query(self, session, pdb):
        """Create a query to look up the rotation matrices.

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get rotation matrices.
        """

        """The query method checks whether data is available in the db
        and if not, recomputes the data. I made the output of the query
        to be empty in order to recompute the data for the unit ids
        without rotation
        """

        rotation = mod.UnitRotations

        test = session.query(rotation).\
                       filter(rotation.unit_id == None)
                       #filter(mod.UnitRotations.unit_id.in_(test_data))

        return test


    def data(self, pdb, **kwargs):
        """Get the rotation matrices for all RNA residues in the given pdb
        that don't already have a rotation matrix.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """

        # The get_units_without_rotation method will get all modified RNA unit_id without rotation
        # based on the list of modified units given in the mod_ref tuple.
        units_without_rotation = self.get_units_without_rotation(pdb)

        # test data:
        # units_without_rotation = ['6N9E|1|1A|5MC|1964', '6N9E|1|1A|5MC|1984', '6N9E|1|1a|5MC|1400', '2L9E|8|A|70U|34', '2L9E|9|A|70U|34', '3T1H|1|X|70U|34', '1FCW|1|D|7MG|46', '1FCW|1|E|7MG|46', '1FIR|1|A|7MG|46', '1GIX|1|B|7MG|46']

        self.logger.info('Calculating rotation matrices for %d nucleotides' % len(units_without_rotation))

        # group_units will be a list of lists of unit ids from one structure at a time
        grouped_units = self.group_units(sorted(units_without_rotation))

        desired_pdb = 'given'    # only update from the list of PDB IDs passed in
        desired_pdb = 'all'      # update from every single PDB ID already in the database

        for unit_id_list in grouped_units:          # unit_id_list is from one structure at a time
            pdb_id = unit_id_list[0].split('|')[0]  # extract PDB id from 0th element of unit_id_list

            if desired_pdb == 'all' or pdb_id in pdb:

                self.logger.info('Loading 3D structure file %s' % pdb_id)
                print('Loading 3D structure file %s' % pdb_id)
                structure = self.structure(pdb_id)      # read .cif file

                for residue in structure.residues():
                    if residue.unit_id() in unit_id_list:
                        if hasattr(residue, 'rotation_matrix'):
                            matrix = residue.rotation_matrix

                            if matrix is not None:
                                print('Adding rotation matrix for %s' % residue.unit_id())
                                self.logger.info('Adding rotation matrix for %s' % residue.unit_id())
                                # Remember to make pdb_id equal to pdb_id and not pdb
                                yield mod.UnitRotations(unit_id=residue.unit_id(),
                                                pdb_id=pdb_id,
                                                cell_0_0=float(matrix[0, 0]),
                                                cell_0_1=float(matrix[0, 1]),
                                                cell_0_2=float(matrix[0, 2]),
                                                cell_1_0=float(matrix[1, 0]),
                                                cell_1_1=float(matrix[1, 1]),
                                                cell_1_2=float(matrix[1, 2]),
                                                cell_2_0=float(matrix[2, 0]),
                                                cell_2_1=float(matrix[2, 1]),
                                                cell_2_2=float(matrix[2, 2]))

