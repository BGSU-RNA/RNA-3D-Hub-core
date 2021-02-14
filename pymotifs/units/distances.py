# coding=utf-8
"""Load distance data for all unit pairs. This will not load distances for
water molecules as those are not generally worth using. The max distance we
will compute distances for is 10.0Å. These values can be changed by altering
the attributes on the Loader class in the module. They cannot be configured in
the configuration file.
"""

import itertools as it
import numpy as np

import pymotifs.core as core
from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader

class Loader(core.SimpleLoader):
    """The actual loader for loading distances.
    """

    # some structures don't have complete nucleotides, otherwise fail to have distances < cutoff
    allow_no_data = True

    max_insert = 5000
    """Number of distances to write at once"""

    dependencies = set([InfoLoader])
    """Stages to depend on"""

    max_distance = 10.0
    """Max distance to use for distances"""

    disallowed = set(['HOH'])
    """Set of components to ignore for distances"""

    def known(self):
        self.logger.info("Querying to find PDBs with distances calculated already")
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id).\
                join(mod.UnitPairsDistances,
                     mod.UnitPairsDistances.unit_id_1 == mod.UnitInfo.unit_id).\
                distinct()
            return set(r.pdb_id for r in query)

    def known_restricted(self,pdbs):
        self.logger.info("Querying to find PDBs with distances calculated already")
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id).\
                join(mod.UnitPairsDistances,
                     mod.UnitPairsDistances.unit_id_1 == mod.UnitInfo.unit_id).\
                filter(mod.UnitInfo.pdb_id.in_(pdbs)).\
                distinct()
            return set(r.pdb_id for r in query)

    def to_process(self, pdbs, **kwargs):
        recalculate = kwargs.get('recalculate', False)
        if recalculate is True:
            return pdbs

        if isinstance(recalculate, (set, list, tuple)):
            try:
                if self.name in recalculate:
                    return pdbs
            except:
                return pdbs

        if bool(self.config[self.name].get('recompute')):
            return True

        if len(pdbs) < 100:
            PDBs_to_process = sorted(set(pdbs) - self.known_restricted(pdbs))
        else:
            PDBs_to_process = sorted(set(pdbs) - self.known())


        if len(PDBs_to_process) > 0:
            return PDBs_to_process
        else:
            raise core.Skip("No new distances to compute")

    def query(self, session, pdb):
        """Create a query to select all entries in unit_pairs_distances for the
        given pdb id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session wrapper to use
        pdb : str
            The pdb id to use.

        Returns
        -------
        query : Query
            A query on unit_pairs_distances for the given pdb.
        """

        with self.session() as session:
            distances = mod.UnitPairsDistances
            return session.query(distances).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == distances.unit_id_1).\
                filter(mod.UnitInfo.pdb_id == pdb)

    def center(self, residue):
        """Get the center to use for the residue. For RNA and AA we use the
        base and backbone centers, otherwise we use an all atom center.

        Parameters
        ----------
        residue : fr3d.data.Component
            The residue to get a center for.

        Returns
        -------
        center : numpy.array
            A numpy array of the center coordinates.

        Note
        ----

        Before November 10, 2020, this program used residue.type to determine
        which kind of center to use, by checking if it equalled 'rna' or 'aa'.
        But those were not actual values of residue.type!
        residue.type apparently comes from the chem_comp.type in the .cif file

        A good way to check is in fr3d-python/fr3d/structures.py, using that

        """

#        self.logger.info('Looking up the center for %s' % residue.unit_id())

        if residue.sequence == 'HOH':
            return None
        if 'base' in residue.centers and len(residue.centers['base']) == 3:      # cover RNA, DNA, and modified nucleotides
            return residue.centers['base']
        if 'aa_fg' in residue.centers and len(residue.centers['aa_fg']) == 3:     # cover amino acids
            return residue.centers['aa_fg']
        if residue.sequence == 'GLY':
            return residue.centers['*']  # same as in fr3d-python\fr3d\data\structures.py

        self.logger.info('No center found for %s, using default coordinates' % residue.unit_id())
#        self.logger.info(np.mean(residue.coordinates(), axis=0))
        self.logger.info(residue.centers['*'])

        return residue.centers['*']  # same as in fr3d-python\fr3d\data\structures.py

    def distance(self, residue1, residue2):
        """Compute the distance between two Components. This will compute the
        distance using the centers from `self.center`. If the Component has no
        center then this will return None.

        Parameters
        ----------
        residue1 : fr3d.data.Component
            The first component
        residue2: fr3d.data.Component
            The second component.

        Returns
        -------
        distance : float, None
            The distance between the residues, or None if there is not a valid
            center for both.
        """

        center1 = self.center(residue1)
        center2 = self.center(residue2)

#        print("distances.py centers of %s and %s" % (residue1.unit_id(),residue2.unit_id()))
#        print(center1)
#        print(center2)

        if center1.size and center2.size:
            return np.linalg.norm(center1 - center2)
        return None

    def is_allowed(self, pair):
        """Check that both entries in the pair are not in self.disallowed.

        Parameters
        ----------
        pair : (fr3d.data.Component, fr3d.data.Component)
            The pair to check.

        Returns
        -------
        allowed : bool
            True if neithe are in `self.disallowed`.
        """

        return pair[0].sequence not in self.disallowed and \
            pair[1].sequence not in self.disallowed

    def data(self, pdb, **kwargs):
        """Compute the distances for all valid pairs of residues in the given
        PDB file. This will not compute distances for things in
        `self.disallowed`.

        Parameters
        ----------
        pdb : str
            The PDB id to use.

        Yields
        ------
        distance : UnitPairDistances
            A UnitPairsDistances entry for the two
        """

        structure = self.structure(pdb)

        # this list of pairs *might* come from fr3d-python/fr3d/structures.py
        # or more likely from fr3d-python/fr3d/data/pairs.py
        # scipy kdtree is used in fr3d-python/fr3d/data/base.py to speed up lookups within a cutoff

        pairs = structure.pairs(distance={'cutoff': self.max_distance})
        pairs = it.ifilter(self.is_allowed, pairs)

        for residue1, residue2 in pairs:
            distance = self.distance(residue1, residue2)
            if distance:
#                self.logger.info("Distance %s %s" % (residue1.unit_id(),residue2.unit_id()))
#                self.logger.info(distance)
                yield mod.UnitPairsDistances(unit_id_1=residue1.unit_id(),
                                             unit_id_2=residue2.unit_id(),
                                             distance=float(distance))

