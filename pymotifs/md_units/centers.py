"""A module to load unit rotation matrices for modified RNA into the database.
"""

from itertools import groupby
from pymotifs import core
from pymotifs import models as mod
from fr3d.modified_parent_mapping import modified_nucleotides


# list the modified nucleotides that should be checked to see if they have a rotation matrix
# if not, one will be computed
# there needs to be a mapping of atoms in each modified nucleotide to the parent atoms,
# in fr3d-python.fr3d.modified_parent_mapping.py
# Apparently this variable needs to be a tuple for the SQL query in get_units_without_rotation
mod_ref = ('DA','DC','DG','DT','A2M', 'CAR', '1MG', '5MC', '5MU', 'H2U', '4SU', 'MIA', 'U34', 'G7M', 'QUO', '1MA', '7MG', 'UDP', 'SRA', 'A23', 'GDP', 'C31', 'U31', 'PGP', 'CCC', '2MA', 'OMC', '2MG', 'M2G', 'OMG', 'YYG', 'YG', '5BU', 'CG1', 'T6A', 'PU', 'PPU', '12A', '2MU', '70U', 'U8U', 'PYY', 'IC', 'IG', 'PSU', 'OMU', 'LCA', 'LCC', 'LCG', 'LKC', 'MAD', 'GAO', 'UAR', '5IC', 'FHU', 'CH', 'CBV', 'IU', 'AVC', '5CG', 'P5P', 'GH3', 'GOM', 'SSU', 'PQ1', 'PUY', 'AET', '0C', '0G', '0U', 'UR3', 'MTU', 'AP7', '31M', '23G', 'MA6', 'MNU', 'A2L', 'C2L', 'G2L', 'U2L', 'RIA', 'A5M', 'M5M', 'UD5', 'S4C', 'FMU', 'PYO', 'SUR', 'P1P', '6IA', 'RTP', '5FU', 'N6G', 'NF2', 'T2T', 'O2C', 'ONE', 'ZAD', 'ZCY', 'ZGU', 'ZHP', 'ZTH', 'A5O', 'RSP', 'RSQ', 'UPV', 'A44', 'G48', 'T39', '2AU', '10C', '4OC', '8AN', 'ZBC', 'ZBU', 'A3P', 'G46', '1SC', 'N5M', 'US5', '1DP', '3TD', '3AU', 'F3O', '5CF', 'A6A', 'A6C', 'A6G', 'A6U', 'GRB', 'RUS', '574', 'DJF', 'C43', 'U36', 'LMS', 'GDO', 'RPC', 'A9Z', 'UBD', 'U37', '1RN', 'A7E', '7AT', '6FC', '6FU', 'URU', '0A', '0U1', '2SG', 'Y5P', '45A', 'M3O', '50L', '50N', '56B', 'U5M', 'CVC', '8AZ', '6F7', '6NW', '6OO', '6OP', 'F2T', '5HM', '73W', '75B', '8OS', '8RJ', 'G4P', '9QV', 'MUM', 'C4J', 'M1Y', 'A7C', 'EQ4', 'U4M', 'UFB', 'UOA', 'UOB', 'B8H', 'B8K', 'B8Q', 'B8T', 'B8W', 'B9B', 'B9H', 'BGH', 'E6G', 'E7G', 'I4U', 'M7A', 'MHG', 'P4U', 'P7G', '4AC', 'B8N', 'E3C', 'RY', 'KGV', 'JMH', 'KAK', 'LHH', 'LV2', 'UY4', 'UY1', 'O2Z', 'U23')
# removed from the list because they are not in modified_parent_mapping:  'TLN', 'ATL', 'N', '3DA', 'C5L', '6MZ', 'I', '5AA', 'F3N', '31H', 'CM0',

# this list was produced directly from modified_parent_mapping, plus DNA nucleotides
mod_ref = ('DA','DC','DG','DT','10C', '125', '126', '127', '12A', '1MA', '1MG', '1SC', '23G', '2AU', '2MA', '2MG', '2MU', '2OM', '3AU', '3MU', '3TD', '4OC', '4SU', '5BU', '5FA', '5FU', '5IC', '5MC', '5MU', '6IA', '70U', '7MG', '8AN', 'A23', 'A2L', 'A2M', 'A3P', 'A44', 'A5M', 'A5O', 'AET', 'AP7', 'AVC', 'C2L', 'C31', 'C43', 'CBV', 'CCC', 'CG1', 'CH', 'CNU', 'CSF', 'FHU', 'G25', 'G2L', 'G46', 'G48', 'G7M', 'GAO', 'GDP', 'GH3', 'GOM', 'GRB', 'GTP', 'H2U', 'IC', 'IG', 'IU', 'KAG', 'LCA', 'M2G', 'M4C', 'M5M', 'MA6', 'MAD', 'MGQ', 'MGV', 'MIA', 'MNU', 'MTU', 'N5M', 'N6G', 'O2G', 'OMC', 'OMG', 'OMU', 'ONE', 'P5P', 'PGP', 'PMT', 'PPU', 'PSU', 'PU', 'PYO', 'QUO', 'RIA', 'RSQ', 'RUS', 'S4C', 'SRA', 'SSU', 'SUR', 'T6A', 'TPG', 'U2L', 'U2P', 'U31', 'U34', 'U36', 'U37', 'U8U', 'UAR', 'UD5', 'UMP', 'UR3', 'URD', 'US5', 'XTS', 'YG', 'YYG', 'ZAD', 'ZBC', 'ZBU', 'ZCY', 'ZGU')

class Loader(core.SimpleLoader):
    """A class to load rotation matrices into the database.
    """

    # dependencies = set([RotationLoader])  # uncomment this if you want to run on new files, which you shouldn't
    allow_no_data = True                    # pipeline won't crash when we pass back no data
    # has_data = False


    # Organize a list of unit ids into a list of lists, one PDB id at a time
    def group_units(self, units_list):

        return [list(i) for j, i in groupby(units_list, lambda a: a.split('|')[0])]


    def get_units_without_centers(self,pdb):
        """ Query the database to find all units with no rotation matrix.
        This will only return units from the tuple mod_ref
        """

        with self.session() as session:

            info = mod.UnitInfo
            centers = mod.UnitCenters

            query = session.query(info.unit_id).\
                            outerjoin(centers, centers.unit_id == info.unit_id).\
                            filter(centers.unit_id == None).\
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

        centers = mod.UnitCenters

        test = session.query(centers).\
                       filter(centers.unit_id == None)
                       #filter(mod.UnitRotations.unit_id.in_(test_data))

        return test


    def data(self, pdb, **kwargs):
        """Get the rotation matrices for all RNA residues in the given pdb
        that don't already have a rotation matrix.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """

        # print out the sequence identifiers for which a mapping to a parent nucleotide is available
        tuple_text = "("
        for k in sorted(modified_nucleotides.keys()):
            tuple_text += "'" + k + "', "
        tuple_text += ")"
        tuple_text = tuple_text.replace(", )",")")
        print(tuple_text)
        self.logger.info("Sequence identifiers with known parent nucleotide: %s" % tuple_text)

        # The get_units_without_rotation method will get all modified RNA unit_id without rotation
        # based on the list of modified units given in the mod_ref tuple.
        units_without_centers = self.get_units_without_centers(pdb)

        grouped_units = self.group_units(sorted(units_without_centers))

        # test_data
        # units_without_centers_list = [['1BZU|1|A|PSU|39'], ['1BZT|1|A|PSU|39'], ['1ASZ|1|R|PSU|655', '1ASZ|1|S|1MG|637', '1ASZ|1|S|5MC|649']]
        # grouped_units = units_without_centers_list

        # self.logger.info('Adding centers data for %d nucleotides' % len(units_without_centers))

        # group_units will be a list of lists of unit ids from one structure at a time
        # grouped_units = self.group_units(sorted(units_without_rotation))

        desired_pdb = 'given'    # only update from the list of PDB IDs passed in
        desired_pdb = 'all'      # update from every single PDB ID already in the database

        self.logger.info("Processing %s 3D structure file(s)" % desired_pdb)

        for unit_id_list in grouped_units:          # unit_id_list is from one structure at a time
            pdb_id = unit_id_list[0].split('|')[0]  # extract PDB id from 0th element of unit_id_list

            if desired_pdb == 'all' or pdb_id in pdb:

                self.logger.info('Loading 3D structure file %s' % pdb_id)
                print('Loading 3D structure file %s' % pdb_id)

                try:
                    structure = self.structure(pdb_id)      # read .cif file

                    foundone = False

                    for residue in structure.residues():
                        if residue.unit_id() in unit_id_list:
                            foundone = True
                            print('Adding centers for %s' % residue.unit_id())
                            self.logger.info('Adding centers for %s' % residue.unit_id())
                            for name in residue.centers.definitions():
                                center = residue.centers[name]
                                if len(center) == 3:
                                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                                          name=name,
                                                          pdb_id=pdb_id,      # use current pdb_id, not pdb the list
                                                          x=float(center[0]),
                                                          y=float(center[1]),
                                                          z=float(center[2]))

                    if not foundone:
                        print("No matches to %s" % unit_id_list)
                        self.logger.info("No matches to %s" % unit_id_list)

                except:
                    print("Could not load %s" % pdb_id)
                    self.logger.info("Could not load %s" % pdb_id)

