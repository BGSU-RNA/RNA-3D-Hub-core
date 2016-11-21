"""Compute a linear ordering for members of an equivalence set.

This will go over all NR classes (equivalence sets) and compute an ordering for
chains that are part of the class. This may skip some chains as we do not
use all chains when computing discrepancies. For example, we skip chains that
have very poor resolution when computing discrepancies. In these cases the
chains will not show up in the final ordering.
"""

import collections as coll

import numpy as np

from sqlalchemy.orm import aliased

from fr3d.ordering.greedyInsertion import orderWithPathLengthFromDistanceMatrix

from pymotifs import core
from pymotifs import models as mod

from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.chains import Loader as NrChainLoader
from pymotifs.nr.classes import Loader as NrClassLoader
from pymotifs.chain_chain.comparision import Loader as SimilarityLoader


class Loader(core.SimpleLoader):
    """The actual Loader to compute and store ordering.
    """

    dependencies = set([NrChainLoader, NrClassLoader, SimilarityLoader])

    trials = 10
    """The number of trials for creating the ordering"""

    def to_process(self, pdbs, **kwargs):
        """Look up all NR classes. This ignores the given PDBs and just creates
        a list of all NR class ids.

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
        classes : list
            A list of all NR class ids.
        """

        latest = None
        if kwargs['manual'].get('nr_release_id', False):
            latest = kwargs['manual']['nr_release_id']
        else:
            data = self.cached(NR_CACHE_NAME)
            if not data:
                raise core.InvalidState("No precomputed grouping to store")
            latest = data['release']

        with self.session() as session:
            query = session.query(mod.NrClasses.nr_class_id).\
                filter_by(nr_release_id=latest)
            return [r.nr_class_id for r in query]

    def query(self, session, class_id):
        """Create a query to find all entries in nr_ordering for the given
        class id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use.

        class_id : int
            The class id

        Returns
        -------
        query : Query
            The query.
        """

        return session.query(mod.NrOrdering).\
            filter_by(nr_class_id=class_id)

    def members(self, class_id):
        """Get all members of the class.

        Parameters
        ----------
        class_id : in
            The id of the the NR class.

        Returns
        -------
        members : list
            A list of tuples (ife_id, nr_chain_id) for all members that are
            part of the class.
        """
        with self.session() as session:
            query = session.query(mod.NrChains.ife_id,
                                  mod.NrChains.nr_chain_id).\
                filter_by(nr_class_id=class_id)
            members = [(r.ife_id, r.nr_chain_id) for r in query]

        if len(members) == 1:
            raise core.Skip("Skip group of size 1")

        if not members:
            raise core.InvalidState("No members in NR class: %i" % class_id)

        return members

    def distances(self, class_id, members):
        """Load all compute distances for members of the NR class. This may not
        load distances for all members, as we do not compute discrepancies for
        all possible chain to chain comparisons. For example, chains with very
        poor resolution are skipped when computing discrepancies.

        Parameters
        ----------
        class_id : int
            The NR class id
        members : list
            A list of members as from `Loader.members`.

        Returns
        -------
        distances : collections.defaultdict
            A dict of dict's that represents the distances. The keys will be
            ife ids, and the values will be the discrepancies between each ife.
        """

        with self.session() as session:
            chains1 = aliased(mod.IfeChains)
            chains2 = aliased(mod.IfeChains)
            nr1 = aliased(mod.NrChains)
            nr2 = aliased(mod.NrChains)
            sim = mod.ChainChainSimilarity
            query = session.query(sim.discrepancy,
                                  chains1.ife_id.label('ife1'),
                                  chains2.ife_id.label('ife2'),
                                  ).\
                join(chains1, chains1.chain_id == sim.chain_id_1).\
                join(chains2, chains2.chain_id == sim.chain_id_2).\
                join(nr1, nr1.ife_id == chains1.ife_id).\
                join(nr2, nr2.ife_id == chains2.ife_id).\
                filter(nr1.nr_class_id == nr2.nr_class_id).\
                filter(nr1.nr_class_id == class_id).\
                order_by(nr1.ife_id, nr2.ife_id)

            distances = coll.defaultdict(lambda: coll.defaultdict(int))
            ifes = set(m[0] for m in members)
            for result in query:
                if result.ife1 not in ifes or result.ife2 not in ifes:
                    continue
                distances[result.ife1][result.ife2] = result.discrepancy

        if not distances:
            raise core.Skip("No distances, skipping class: %i" % class_id)

        return distances

    def ordered(self, members, distances):
        """Compute an ordering for the members of an equivalence set given the
        distance matrix.

        Parameters
        ----------
        members : list
            A list of members as from `Loader.members`.
        distances : collections.defaultdict
            A defaultdict of distances as from `Loader.distances`.

        Returns
        -------
        ordered_members : list
            A list of the given members in an order as specified by the
            discrepancies. This will not contain all members of the given list
            as those things without distances will be skipped.
        """

        dist = np.zeros((len(members), len(members)))
        for index1, member1 in enumerate(members):
            curr = distances.get(member1, {})
            for index2, member2 in enumerate(members):
                val = curr.get(member2, None)
                if member2 not in curr:
                    val = None
                dist[index1, index2] = val

        ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist, self.trials)
        return [members[index] for index in ordering]

    def data(self, class_id, **kwargs):
        """Compute the ordering rows to store. This will lookup the distances
        and members as needed. It is possible for this to raise a
        `pymotifs.core.Skip` if there are no distances stored.

        Parameters
        ----------
        class_id : int
            The class id to compute distances for.

        Raises
        ------
        core.Skip
            If there is no distances.

        Returns
        -------
        data : list
            A list of NrOrdering objects to store.
        """

        members = self.members(class_id)
        distances = self.distances(class_id, members)
        ordered = self.ordered(members, distances)
        data = []
        for index, (ife_id, chain_id) in enumerate(ordered):
            data.append(mod.NrOrdering(
                nr_chain_id=chain_id,
                nr_class_id=class_id,
                index=index,
            ))
        return data
