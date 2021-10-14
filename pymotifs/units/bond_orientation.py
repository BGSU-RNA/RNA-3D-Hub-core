"""A module to compute and store glycosidic bond orientations in the database.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader

import numpy as np
import math
from fr3d.modified_parent_mapping import modified_nucleotides


class Loader(core.SimpleLoader):
    """A class to compute glycosidic bond orientations and load them into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def query(self, session, pdb):
        """Create a query to look up glycosidic bond orientations

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get bond orientations.
        """

        return session.query(mod.UnitOrientations).\
            filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        """
        Load the structure, compute the orientations, save.

        :pdb: The pdb to process
        :yields: Yields angle chi_degree and an annotation
        """

        # load the 3D structure file
        structure = self.structure(pdb)

        nts = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides

        for nt in nts:
            chi = None
            orientation = ""

            N1N9 = []
            C2C4 = []

            if nt.sequence in ['A','G','DA','DG']:
                N1N9 = nt.centers["N9"]
                C2C4 = nt.centers["C4"]
            elif nt.sequence in ['C','U','DC','DT']:
                N1N9 = nt.centers["N1"]
                C2C4 = nt.centers["C2"]
            elif nt.sequence in modified_nucleotides:
                parent = modified_nucleotides[nt.sequence]["standard"]
                if parent in ['A','G','DA','DG']:
                    N1N9 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N9"]]
                    C2C4 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["C4"]]
                elif parent in ['C','U','DC','DT']:
                    N1N9 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N1"]]
                    C2C4 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["C2"]]
                else:
                    print("%s has no identified parent" % nt.unit_id())
                    self.logger.info("%s has no identified parent" % nt.unit_id())

            if len(N1N9) == 3 and len(C2C4) == 3:

                try:
                    O4P_C1P = nt.centers["C1'"] - nt.centers["O4'"]
                    N1N9_C1P  = nt.centers["C1'"] - N1N9
                    C2C4_N1N9 = N1N9 - C2C4

                    # chi angle definition:
                    # O4*_C1*_N1_C2 (for pyrimidines)
                    # O4*_C1*_N9_C4 (for purines):

                    perp_to_sugar       = np.cross(O4P_C1P,N1N9_C1P)
                    norm_perp_to_sugar  = np.linalg.norm(perp_to_sugar)
                    if norm_perp_to_sugar != 0:
                        perp_to_sugar   = perp_to_sugar/norm_perp_to_sugar

                    perp_to_base        = np.cross(-N1N9_C1P,C2C4_N1N9)
                    norm_perp_to_base   = np.linalg.norm(perp_to_base)
                    perp_to_base        = perp_to_base/norm_perp_to_base

                    cross_cross_chi = np.cross(perp_to_base,perp_to_sugar)

                    # Take the dot product of the vectors perp_to_base &
                    # perp_to_sugar to get cos(chi).
                    # Take norm(cross product) to get sin(chi).

                    cos_chi = np.dot(perp_to_sugar,perp_to_base)
                    if np.dot(cross_cross_chi,N1N9_C1P) > 0:
                        sin_chi = np.linalg.norm(cross_cross_chi)
                    else:
                        sin_chi = -np.linalg.norm(cross_cross_chi)

                    # sign of chi_degree matches Bevilacqua 2011 paper on syn and anti
                    # this definition matches the IUPAC definition from http://www.chem.qmul.ac.uk/iupac/misc/pnuc2.html#230

                    chi  = 180*math.atan2(sin_chi,cos_chi)/math.pi # glycosidic bond angle

                    # Giving nomenclature according to chi values: anti (most common), or syn
                    if chi > -90 and chi < -45:
                        orientation = 'Intermediate Syn'
                    elif chi >= -45 and chi < 90:
                        orientation = 'Syn'
                    else:
                        orientation = 'Anti'

                    if chi and orientation:
                        yield mod.UnitOrientations(unit_id=nt.unit_id(),
                                                pdb_id=pdb,
                                                orientation=orientation,
                                                chi_degree=chi)

                except:
                    self.logger.info("%s glycosidic bond orientation calculation failed" % nt.unit_id())