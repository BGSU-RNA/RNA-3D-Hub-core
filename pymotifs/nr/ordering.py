"""Compute a linear ordering for members of an equivalence set.

This will go over all NR classes (equivalence sets) and compute an ordering for
chains that are part of the class. This may skip some chains as we do not
use all chains when computing discrepancies. For example, we skip chains that
have very poor resolution when computing discrepancies. In these cases the
chains will not show up in the final ordering.
"""

import collections as coll
from collections import defaultdict

import numpy as np
import time

from sqlalchemy import func, select
from sqlalchemy.orm import aliased

from fr3d.ordering.greedyInsertion import orderWithPathLengthFromDistanceMatrix
from orderBySimilarity import imputeNANValues
from orderBySimilarity import treePenalizedPathLength

from pprint import pprint

from pymotifs import core
from pymotifs import models as mod

from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.chains import Loader as NrChainLoader
from pymotifs.nr.classes import Loader as NrClassLoader
from pymotifs.nr.cqs import NrQualityLoader
from pymotifs.chain_chain.comparision import Loader as SimilarityLoader


class Loader(core.SimpleLoader):
    """The actual Loader to compute and store ordering.

    Attributes
    ----------
    trials : int, 10
        The number to runs to use when sorting the members of each group.
    """

    trials = 100
    dependencies = set([NrChainLoader, NrClassLoader, NrQualityLoader, SimilarityLoader])

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

        # determine desired release or retrieve most recent release
        latest = None
        if kwargs.get('manual', {}).get('nr_release_id', False):
            latest = kwargs['manual']['nr_release_id']
        else:
            data = self.cached(NR_CACHE_NAME)
            if not data:
                raise core.InvalidState("No precomputed grouping to store")
            latest = data['release']

        self.logger.info("to_process: latest: %s" % latest)

        # find all nr_class_id values associated with this release
        if 0 > 1:
            with self.session() as session:
                query = session.query(mod.NrClasses.nr_class_id).\
                    filter_by(nr_release_id=latest)
                return [(latest, r.nr_class_id) for r in query]

        # fill in all orderings from release 2.* and 3.*
        with self.session() as session:
            # get all triples of nr_release_id and nr_class_id and equivalence class name
            query = session.query(mod.NrClasses.nr_release_id,mod.NrClasses.nr_class_id,mod.NrClasses.name)

            # retrieve all nr_class_id values after releases 0.* and 1.*
            all_triples = []
            for r in query:
                if not r.nr_release_id[0:2] == "0." and not r.nr_release_id[0:2] == "1.":
                    all_triples.append((r.nr_release_id,r.nr_class_id,r.name))

            # sort by equivalence class name and then by release number
            all_triples = sorted(all_triples, key = lambda r: (r[2], r[1]))

            self.logger.info("to_process: sorted triples", all_triples)

            # get all pairs of nr_class_id and ife_id
            query = session.query(mod.NrOrderingTest.nr_class_name)
            ordered_nr_class_name = [r.nr_class_name for r in query]

            # loop through triples, save first class_id for each name, save pair if not already ordered
            one_class_id_per_name = {}
            pairs_to_process = []
            for (nr_release,nr_class,nr_name) in all_triples:
                if not nr_name in one_class_id_per_name:
                    one_class_id_per_name[nr_name] = (nr_release,nr_class)
                    if not nr_name in ordered_nr_class_name:
                        self.logger.info("to_process: need to order release %s class %s" % (nr_release,nr_name))
                        pairs_to_process.append((nr_release,nr_class))
                    else:
                        self.logger.info("to_process: already ordered release %s class %s" % (nr_release,nr_name))

            return pairs_to_process

#        return [(latest, r.nr_class_id) for r in query]


            # the pairs returned here will be passed to the data method one by one
            # make sure to only pass each nr_class_id one time, not once for each release


    #def is_missing(self, entry, **kwargs):
    #    """Placeholder to see how to properly ID the classes that need
    #    attention.
    #    """
    #    pass

    def has_data(self, entry, **kwargs):
        # for each equivalence class id, which is an integer,
        # check if there is already an ordering with enough data in it
        # release rel is not used
        rel, class_id = entry

        with self.session() as session:
            # find full equivalence class names matching the class_id in entry
            # NrClasses is the table nr_classes
            # .name is something like NR_1.5_01181.1
            query = session.query(mod.NrClasses.name).\
                filter_by(nr_class_id=class_id)

            # the list below should have length one, usually does
            # only the 0 element will be used
            class_name = [r.name for r in query]

            # NrChains is the table nr_chains
            query = session.query(mod.NrChains.ife_id).\
                filter_by(nr_class_id=class_id)

            # list of the ife_id's associated with this equivalence class
            results = [r.ife_id for r in query]

            # retrieve all ordering information for the first class_name
            query = session.query(mod.NrOrderingTest.ife_id).filter_by(nr_class_name=class_name[0])

            count = query.count()

            self.logger.info("has_data: %s previously-ordered IFEs in class %s (id: %s, %s members): %s" %
                             (count, class_name[0], class_id, len(results), str(results)))

            if not count == len(results):
                self.logger.info("has_data: Number previously ordered is different from number of members.")

            # if there is one piece of ordering data for each IFE in this class, we have data
            # if count == len(results):  # a small number of classes have different numbers; don't re-order
            if count > 0:
                return True

        return False

    def query(self, session, pair):
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

        _, class_id = pair

        self.logger.info("query: class_id: %s" % class_id)

        return session.query(mod.NrOrdering).\
            filter_by(nr_class_id=class_id)

    def get_original_info(self, class_name):
        """Obtain the first nr_release_id and the corresponding
        nr_class_id value for the input class_name.

        Parameters
        ----------
        class_name : in
            The name of the the NR class.

        Returns
        -------
        orig_release_id : string
            The first release_id in which the NR class appears.
        """

        with self.session() as session:
            ncl = aliased(mod.NrClasses)
            nre = aliased(mod.NrReleases)

            query = session.query(nre.nr_release_id,
                                  ncl.nr_class_id).\
                join(ncl, ncl.nr_release_id == nre.nr_release_id).\
                filter(ncl.name == class_name).\
                order_by(nre.index).\
                limit(1)

            for r in query:
                orig_release_id = r.nr_release_id
                orig_class_id = r.nr_class_id

            return (orig_release_id, orig_class_id)

    def members_revised(self, class_id, release_id):
        """Get all members of the class.

        Parameters
        ----------
        class_id : in
            The first class_id value for the NR class.

        release_id : in
            The first representative sets release that contains the class.

        Returns
        -------
        members : list
            A list of tuples (ife_id, nr_chain_id) for all
            members of the class.
        """

        self.logger.info("members_revised:  class_id: %s" % class_id)

        with self.session() as session:
            nch = aliased(mod.NrChains)

            query = session.query(nch.ife_id, nch.nr_chain_id).\
                filter(nch.nr_class_id == class_id)

            members = [(r.ife_id, r.nr_chain_id) for r in query]

#        if len(members) == 1:
#            raise core.Skip("Skip group of size 1")

        if not members:
            raise core.InvalidState("No members in NR class: %i" % class_id)

        return members

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
            ife_id is like 2A43|1|A and nr_chain_id is like 11890928
        """

        self.logger.info("members: class_id: %s" % class_id)

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

    def distances_revised(self, release_id, class_id, members):
        """Load all available computed distances between members of the NR class.
        Note that all possible chain-to-chain comparisons are not computed;
        for example, chains with very poor resolution are skipped during the
        discrepancy calculations.

        Parameters
        ----------
        class_id : in
            The first class_id value for the NR class.

        release_id : in
            The first representative sets release that contains the class.

        members : list
            A list of members as from `Loader.members`.  (Actually generated
            in members_revised)

        Returns
        -------
        distances_revised : collections.defaultdict
            A dict-of-dicts that represents the distances. The keys will be
            ife_ids, and the values will be the discrepancies between each pair
            of IFEs.
        """

        self.logger.info("distances_revised: class_id (%s) has %s members including %s" % (class_id, len(members), members[0]))

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
                filter(nr1.nr_release_id == nr2.nr_release_id).\
                filter(nr1.nr_release_id == release_id).\
                order_by(nr1.ife_id, nr2.ife_id)

            distances_revised = coll.defaultdict(lambda: coll.defaultdict(int))

            ifes = set(m[0] for m in members)

            for result in query:
                if result.ife1 not in ifes or result.ife2 not in ifes:
                    continue
                distances_revised[result.ife1][result.ife2] = result.discrepancy

        if not distances_revised:
            raise core.Skip("No distances, skipping class: %i" % class_id)

        if set(distances_revised.keys()) != ifes:
            missing = ', '.join(ifes - set(distances_revised.keys()))
            self.logger.warning("Did not load distances for all pairs in: %i."
                                " Missing %s", class_id, missing)

        return distances_revised

    def distances(self, nr_release_id, class_id, members):
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

        self.logger.info("distances: class_id (%s) has %s members" % (class_id, len(members)))

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
                filter(nr1.nr_release_id == nr2.nr_release_id).\
                filter(nr1.nr_release_id == nr_release_id).\
                order_by(nr1.ife_id, nr2.ife_id)

            distances = coll.defaultdict(lambda: coll.defaultdict(int))
            ifes = set(m[0] for m in members)
            for result in query:
                if result.ife1 not in ifes or result.ife2 not in ifes:
                    continue
                distances[result.ife1][result.ife2] = result.discrepancy

        if not distances:
            raise core.Skip("No distances, skipping class: %i" % class_id)

        if set(distances.keys()) != ifes:
            missing = ', '.join(ifes - set(distances.keys()))
            self.logger.warning("Did not load distances for all pairs in: %i."
                                " Missing %s", class_id, missing)

        return distances

    def ordered_revised(self, members, distances):
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
            A sparse list of the given members in an order as specified by the
            discrepancies. (Sparse because any members without distances will
            be skipped.)
        """

        self.logger.info("ordered_revised: %s members" % len(members))

        dist = np.zeros((len(members), len(members)))

        for index1, member1 in enumerate(members):
            curr = distances.get(member1[0], {})
            for index2, member2 in enumerate(members):
                val = curr.get(member2[0], None)
                if member2[0] not in curr:
                    val = None
                dist[index1, index2] = val
                self.logger.debug("ordered: dist[%s, %s] = %s" % (index1, index2, val))

        newDist = imputeNANValues(dist)
        print(newDist)
        ordering = treePenalizedPathLength(newDist,max(self.trials,len(members)))

        return [members[index] for index in ordering]

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

        self.logger.info("ordered: %s members" % len(members))

        dist = np.zeros((len(members), len(members)))
        for index1, member1 in enumerate(members):
            curr = distances.get(member1[0], {})
            for index2, member2 in enumerate(members):
                val = curr.get(member2[0], None)
                if member2[0] not in curr:
                    val = None
                dist[index1, index2] = val
                self.logger.debug("ordered: dist[%s, %s] = %s" % (index1, index2, val))

        ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist,
                                                               self.trials,
                                                               scanForNan=True)
        return [members[index] for index in ordering]

    def mark_processed(self, pair, **kwargs):
        return super(Loader, self).mark_processed(pair[1], **kwargs)

    def get_nrclassname(self, class_id):
        with self.session() as session:
            query = session.query(mod.NrClasses.name).\
                filter_by(nr_class_id=class_id)

            return [r.name for r in query]

    def data(self, pair, **kwargs):
        """Compute the ordering rows to store. This will look up the distances
        and members as needed. It is possible for this to raise a
        `pymotifs.core.Skip` if there are no distances stored.

        Parameters
        ----------
        pair consists of nr_release_id and class_id

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

        nr_release_id, class_id = pair

        # query database to get class name like NR_1.5_09579.1
        get_nr_class_name = self.get_nrclassname(class_id)
        nr_class_name = get_nr_class_name[0]

        self.logger.info("data: INPUT: nr_release_id: %s, class_id: %s, nr_class_name: %s" % (nr_release_id, class_id, nr_class_name))

        # look up release_id and class_id when this class was first created, based on its name
        # that way we can just focus on one class_id instead of keeping all of them in play
        orig_release_id, orig_class_id = self.get_original_info(nr_class_name)
        self.logger.info("data: USING: orig_release_id %s and orig_class_id %s for nr_class_name %s for class_id %s"
                         % (orig_release_id, orig_class_id, nr_class_name, class_id))

        # members = self.members(class_id)
        # self.logger.info("data: members: %s (class_id %s)" % (str(members), class_id))

        members_revised = self.members_revised(orig_class_id, orig_release_id)
        self.logger.info("data: members_revised: %s (class_id %s)" % (str(members_revised), orig_class_id))

        if len(members_revised) <= 2:
            # no need to try to find an ordering, all possible orderings are equivalent
            ordered_revised = members_revised
        elif len(members_revised) <= 300:
            # look up distances using database for smallish groups
            # on 10/15/2019, database lookup of a group with 299 members took 1.04 seconds
            # flat file reading of groups up to 450 members took under 2 seconds
            # 300 is a good cutoff between the two
            # before reading the flat file, it could take hours to look up discrepancies from the database for large groups

            starttime = time.clock()
            #distances = self.distances(nr_release_id, class_id, members)
            distances_revised = self.distances_revised(orig_release_id, orig_class_id, members_revised)

            self.logger.info("data: time to read a group of size %d from database was %8.4f seconds" % (len(members_revised),time.clock()-starttime))

            #self.logger.info("data: distances: %s (class_id %s)" % (str(distances), class_id))
            #self.logger.info("data: distances_revised: %s (class_id %s)" % (str(distances_revised), orig_class_id))

            #ordered = self.ordered(members, distances)
            ordered_revised = self.ordered_revised(members_revised, distances_revised)

        else:
            self.logger.info("large group %s:  reading flat file of discrepancies" % nr_class_name)
            print("large group:  reading flat file of discrepancies")
            discDict = defaultdict(lambda: None)
            # put discrepancy information into a dictionary
            starttime = time.clock()
            infileNameWithPath = "/tmp/discrepancy2.txt"
            with open(infileNameWithPath,"r") as infile:
                for line in infile:
                    fields = line.replace("\n","").split("\t")
                    ife1 = fields[0]
                    ife2 = fields[1]
                    discrepancy = float(fields[2])
                    discDict[(ife1,ife2)] = discrepancy
                    discDict[(ife2,ife1)] = discrepancy

            self.logger.info("large group:  read flat file of discrepancies")
            self.logger.info("data: time to read a group of size %d from flat file was %8.4f seconds" % (len(members_revised),time.clock()-starttime))

            dist = np.zeros((len(members_revised), len(members_revised)))

            for index1, member1 in enumerate(members_revised):
                for index2, member2 in enumerate(members_revised):
                    val = discDict[(member1[0],member2[0])]
                    dist[index1, index2] = val

            newDist = imputeNANValues(dist)
            print(newDist)
            ordering = treePenalizedPathLength(newDist,max(self.trials,len(members_revised)))

            ordered_revised = [members_revised[index] for index in ordering]
            self.logger.info("large group:  produced a new ordering")


        #self.logger.info("data: ordered: %s (class_id %s)" % (str(ordered), class_id))
        #self.logger.info("data: ordered_revised: %s (class_id %s)" % (str(ordered_revised), orig_class_id))

        #data = []
        #for index, (ife_id, chain_id) in enumerate(ordered):
        #    self.logger.info("data: data: class_id %s / index %s / chain_id %s" % (class_id, index, chain_id))
        #    data.append(mod.NrOrdering(
        #        nr_chain_id=chain_id,
        #        nr_class_id=class_id,
        #        index=index,
        #    ))

        data_revised = []
        for index, (ife_id, nr_chain_id) in enumerate(ordered_revised):
            self.logger.info("data: data_revised: nr_class_id %s / index %s / nr_chain_id %s / nr_class_name %s / ife_id %s " % (orig_class_id, index, nr_chain_id, nr_class_name, ife_id))
            data_revised.append(mod.NrOrderingTest(
                nr_class_id=orig_class_id,
                nr_class_name=nr_class_name,
                nr_chain_id=nr_chain_id,
                ife_id=ife_id,
                class_order=index,
            ))

        #self.logger.info("data: output data: %s (class_id %s)" % (repr(data), class_id))
        #self.logger.info("data: output data_revised: %s (class_id %s)" % (repr(data_revised), class_id))

        return data_revised
        #return data
