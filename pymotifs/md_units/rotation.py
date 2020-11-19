"""A module to load unit rotation matrices for modified RNA into the database.
"""

from pymotifs import core
from pymotifs import models as mod

mod_ref = ('A2M', 'CAR', '1MG', '5MC', '5MU', 'H2U', '4SU', 'MIA', 'U34', 'G7M', 'QUO', '1MA', '7MG', 'UDP', 'SRA', 'A23', 'GDP', 'C31', 'U31', 'PGP', 'CCC', 'N', '2MA', 'OMC', '2MG', 'M2G', 'OMG', 'YYG', 'YG', '5BU', 'CG1', 'I', 'T6A', 'PU', 'PPU', '12A', '2MU', '70U', 'U8U', 'PYY', 'IC', 'IG', 'PSU', 'OMU', 'ATL', 'LCA', 'LCC', 'LCG', 'LKC', 'TLN', 'MAD', 'GAO', 'UAR', '5IC', 'FHU', 'CH', 'CBV', 'IU', 'AVC', '5CG', 'P5P', 'GH3', 'GOM', 'SSU', 'PQ1', 'PUY', 'AET', '0C', '0G', '0U', 'UR3', '3DA', 'MTU', 'AP7', '5AA', 'F3N', '31H', '31M', '23G', 'MA6', 'MNU', 'A2L', 'C2L', 'G2L', 'U2L', 'RIA', 'A5M', 'M5M', 'UD5', 'S4C', 'FMU', 'PYO', 'SUR', 'P1P', '6IA', 'RTP', '5FU', 'N6G', 'NF2', 'T2T', 'O2C', 'ONE', 'ZAD', 'ZCY', 'ZGU', 'ZHP', 'ZTH', '6MZ', 'CM0', 'A5O', 'RSP', 'RSQ', 'UPV', 'A44', 'C5L', 'G48', 'T39', '2AU', '10C', '4OC', '8AN', 'ZBC', 'ZBU', 'A3P', 'G46', '1SC', 'N5M', 'US5', '1DP', '3TD', '3AU', 'F3O', '5CF', 'A6A', 'A6C', 'A6G', 'A6U', 'GRB', 'RUS', '574', 'DJF', 'C43', 'U36', 'LMS', 'GDO', 'RPC', 'A9Z', 'UBD', 'U37', '1RN', 'A7E', '7AT', '6FC', '6FU', 'URU', '0A', '0U1', '2SG', 'Y5P', '45A', 'M3O', '50L', '50N', '56B', 'U5M', 'CVC', '8AZ', '6F7', '6NW', '6OO', '6OP', 'F2T', '5HM', '73W', '75B', '8OS', '8RJ', 'G4P', '9QV', 'MUM', 'C4J', 'M1Y', 'A7C', 'EQ4', 'U4M', 'UFB', 'UOA', 'UOB', 'B8H', 'B8K', 'B8Q', 'B8T', 'B8W', 'B9B', 'B9H', 'BGH', 'E6G', 'E7G', 'I4U', 'M7A', 'MHG', 'P4U', 'P7G', '4AC', 'B8N', 'E3C', 'RY', 'KGV', 'JMH', 'KAK', 'LHH', 'LV2', 'UY4', 'UY1', 'O2Z', 'U23')

#from pymotifs.units.info import Loader as InfoLoader

class Loader(core.SimpleLoader):
    """A class to load rotation matrices into the database.
    """

    #dependencies = set([InfoLoader])
    #allow_no_data = True
    # has_data = False

    def get_modified_nts(self):
        
        with self.session() as session:

            info = mod.UnitInfo
            rotation = mod.UnitRotations
            
            query = session.query(info.unit_id).\
                            outerjoin(rotation, rotation.unit_id == info.unit_id).\
                            filter(rotation.unit_id == None).\
                            filter(info.unit_type_id == None).\
                            filter(info.unit.in_(mod_ref))

        return [r.unit_id for r in query]
        
    
    def query(self, session, pdb):
        """Create a query to lookup the rotation matrices.

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
        """Get the rotation matrices for all RNA residues in the given pdb.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """

        test_data = ['5JTE|1|AW|7MG|46', '5JTE|1|AW|3AU|47', '5JTE|1|AW|5MU|54']

        # The get_modified_nts method will get all modified RNA unit_id without rotation
        # based on the list of modified units given in the mod_ref tuple.
        test_mod = self.get_modified_nts()
        print test_mod
        print len(test_mod)

