"""Compute a summary of the number for each type of interaction for each unit.
This stage will compute the number of each type of basepair, base stack, and
base phosphate (with special cases for long range) that each unit in a
structure is the first entry of. This is useful in fr3d for searching as well
as for general querying and summarizing.
"""

import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import LONG_RANGE
from pymotifs.interactions.pairwise import IGNORE

from pymotifs.interactions.pairwise import Loader as InterLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader


class Loader(core.SimpleLoader):
    dependencies = set([InterLoader, UnitLoader, PdbLoader])

    ignore_bp = IGNORE
    """A list of basepair families to ignore the counts of."""

    @property
    def table(self):
        return mod.UnitInteractionSummary

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions.pdb_id).\
                distinct()
            known = set(r.pdb_id for r in query)

        return sorted(set(pdbs).intersection(known))

    def query(self, session, pdb):
        """Build a query to find all summary entries for the given PDB.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session to use
        pdb : str
            The pdb id to use

        Returns
        -------
        query : Query
            The query to use.
        """

        return session.query(mod.UnitInteractionSummary).\
            filter_by(pdb_id=pdb)

    def increment_bp(self, current, bp, crossing):
        """Increment the count of the current bp. If the base pair is in
        `Loader.ignore_bp` we return the current counts.

        Parameters
        ----------
        current : dict
            The current dictonary of counts
        bp : str
            The type of basepair to increment
        crossing : int
            The crossing number to increment

        Returns
        -------
        counts : dict
            The updated counts.
        """

        if bp in self.ignore_bp:
            return current
        return self.increment(current, 'bps', bp, crossing)

    def increment_stacks(self, current, stack, crossing):
        """Increment the counts for the given stack annotation.

        Parameters
        ----------
        current : dict
            The current counts
        stack : str
            The stacking annotation
        crossing : int
            The crossing number

        Returns
        -------
        counts : dict
            The updated number of counts.
        """
        return self.increment(current, 'stacks', stack, crossing)

    def increment_bphs(self, current, unit1, unit2, bph, crossing):
        """Increment the counts for the given base phosphate annotation. We do
        not count self 0BPh interactions as those are very common.

        Parameters
        ----------
        current : dict
            The current counts
        bph : str
            The base phosphaate annotation
        crossing : int
            The crossing number

        Returns
        -------
        counts : dict
            The updated number of counts.
        """
        if unit1 != unit2 and bph != '0BPh':
            return self.increment(current, 'bphs', bph, crossing)
        return current

    def increment(self, current, family, name, crossing):
        """Increment the the counts of the given annotation for the given name.
        This will increment totals as well as the long range counts if the
        given annotation is long range. We do not increment counts if the
        interaction is near which we can tell by it starting with 'n'.

        Parameters
        ----------
        current : dict
            The current counts
        family : str
            The family of interaction, like 'bp', or 'bph', etc
        name : str
            The annotation to increment
        crossing : int
            The corssing number.
        """

        if name and name[0] != 'n':
            current[name] += 1
            current['total'] += 1
            current[family] += 1
            lr_inc = int(crossing > LONG_RANGE)
            current['lr_' + name] += lr_inc
            current['lr_total'] += lr_inc
            current['lr_' + family] += lr_inc
        return current

    def data(self, pdb_id, **kwargs):
        """Compute the summary for all units in the given pdb. This will look
        up all RNA bases in the given structure and compute a summary of the
        number of interactions for each unit.

        Parameters
        ----------
        pdb_id : str
            The pdb id.

        Returns
        -------
        summaries : list
            A list of dictonaries as from 'increment'.
        """

        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id.label('unit_id_1'),
                                  mod.UnitInfo.model,
                                  mod.UnitInfo.chain,
                                  mod.UnitInfo.pdb_id,
                                  mod.UnitPairsInteractions.unit_id_2,
                                  mod.UnitPairsInteractions.f_lwbp,
                                  mod.UnitPairsInteractions.f_bphs,
                                  mod.UnitPairsInteractions.f_stacks,
                                  mod.UnitPairsInteractions.f_crossing,
                                  ).\
                outerjoin(mod.UnitPairsInteractions,
                          mod.UnitInfo.unit_id == mod.UnitPairsInteractions.unit_id_1).\
                filter(mod.UnitInfo.pdb_id == pdb_id).\
                filter(mod.UnitInfo.unit_type_id == 'rna')

            data = coll.defaultdict(lambda: coll.defaultdict(int))
            for result in query:
                current = data[result.unit_id_1]
                current['unit_id'] = result.unit_id_1
                current['pdb_id'] = result.pdb_id
                current['model'] = result.model
                current['chain'] = result.chain
                crossing = result.f_crossing
                self.increment_bp(current, result.f_lwbp, crossing)
                self.increment_stacks(current, result.f_stacks, crossing)
                self.increment_bphs(current, result.unit_id_1,
                                    result.unit_id_2, result.f_bphs, crossing)
                data[current['unit_id']] = current

            return [(dict(v)) for v in data.values()]
