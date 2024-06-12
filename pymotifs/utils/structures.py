"""Some queries that are useful in accessing stuff in the database.
"""

import itertools as it

from sqlalchemy.orm import aliased

from pymotifs import core

from pymotifs import utils as ut
from pymotifs import models as mod

SYNTHEIC = (32630, 'synthetic construct')

LONG_RANGE_CUTOFF = 4

CURRENT_REP_QUERY = '''
select rep_id
from nr_current_representative
where
    pdb_id = :val
'''

CURRENT_MEMBERS_QUERY = '''
select N2.pdb_id
from nr_current_representative as N1
join nr_current_representative as N2
on
    N2.rep_id = N1.rep_id
where
    N1.pdb_id = :val
;
'''


class UnknownTaxonomyException(Exception):
    """This is raised when we are looking for taxonomy ids but cannot find them
    all.
    """
    pass


class Base(object):
    """A base class for the query classes.
    """
    def __init__(self, maker):
        self.session = core.Session(maker)


class NR(Base):
    def members(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_MEMBERS_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep

    def representative(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_REP_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep


class Structure(Base):
    def rna_chains(self, pdb, return_id=False, strict=False, extended=True):
        """This will get all chains labeled as RNA for a given structure or
        structures. This has a strict mode which can filter out chains which
        are not standard RNA, however, this may also filter out chains where
        are RNA incorrectly. For example, things with modified bases listed in
        their sequence. The strict mode is needed when dealing with things like
        3CPW|A which is labeled as RNA but is actually a protein chain.

        :pdb: The pdb id or a list of pdb ids to query.
        :return_id: A boolean to indicate if we should return a tuple that is
        chain_id, chain_name or just chain_name.
        :strict: A flag that will exclude any chains that are not composed of
        only ACGUN.
	:extended: A flag that sets the allowed macromolecule types for
	processing.  False = "Polyribonucleotide (RNA)" only; True allows
	additional values.
        :returns: A list of the names or a tuple of the ids and names.
        """

        # in November 2020, the type for RNA comes back as polyribonucleotide
        # note:  this is for entity_macromolecule_type, see below
        # note:  it would be really helpful to find a less fragile way to do this!
        # also set in mapping.py
        macromolecule_types = set(['Polyribonucleotide (RNA)','polyribonucleotide'])
#        if extended:
        macromolecule_types.add('DNA/RNA Hybrid')
        macromolecule_types.add('NA-hybrid')
        macromolecule_types.add('polydeoxyribonucleotide/polyribonucleotide hybrid')



        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name,
                     mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(macromolecule_types))

            if isinstance(pdb, basestring):
                query = query.filter_by(pdb_id=pdb)
            else:
                query = query.filter(mod.ChainInfo.pdb_id.in_(pdb))

            # look specifically for A, C, G, U, N but not modified nts
            if strict:
                func = mod.ChainInfo.sequence.op('regexp')
                query = query.filter(func('^[ACGUN]$'))

            data = []
            for result in query:
                entry = result.chain_name
                if return_id:
                    entry = (result.chain_name, result.chain_id)
                data.append(entry)

            return data

    def longest_chain(self, pdb, model=1):
        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter_by(pdb_id=pdb, model=model).\
                filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(mod.UnitInfo.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def interactions(self, pdb, chain):
        c1 = aliased(mod.UnitInfo)
        c2 = aliased(mod.UnitInfo)
        interactions = []
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions).\
                join(c1, c1.id == mod.UnitPairsInteractions.unit_id_1).\
                join(c2, c2.id == mod.UnitPairsInteractions.unit_id_2).\
                filter(mod.UnitPairsInteractions.pdb_id == pdb).\
                filter(c1.chain == c2.chain, c1.chain == chain)

            for result in query:
                data = ut.row2dict(result)
                data['id'] = int(data['unit_pairs_interactions_id'])
                interactions.append(data)

        return interactions

    def stacks(self, pdb, chain, count=None):
        """Get the list or count of all stacking interactions.

        :pdb: The pdb file to query.
        :chain: The chain to get stacks for.
        :count
        """
        pass

    def loops(self, pdb):
        loops = []
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.LoopsAll, mod.LoopPositions.loop_id == mod.LoopsAll.loop_id).\
                filter(mod.LoopsAll.pdb_id == pdb).\
                order_by(mod.LoopsAll.loop_id, mod.LoopPositions.position)

            grouped = it.groupby(it.imap(ut.row2dict, query),
                                 lambda a: a['loop_id'])

            for loop_id, positions in grouped:
                positions = list(positions)
                endpoints = []
                for pos in positions:
                    if pos.border:
                        endpoints.append(pos['unit_id'])

                loops.append({
                    'id': loop_id,
                    'nts': [pos['unit_id'] for pos in positions],
                    'endpoints': list(ut.grouper(2, endpoints))
                })

        return loops

    def source(self, pdb, chain, simplify=False):
        """This is a method to extract all species level taxonomy ids for a
        given chain. A chain can be composed of more than one, thus a list.
        However, in some cases there is no known taxnomy id, in this case an
        empty list is returned.

        :pdb: The pdb to query.
        :chain: The chain to query.
        :returns: A list of taxonomy ids the chain has.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.taxonomy_id).\
                filter_by(pdb_id=pdb).\
                filter_by(chain_name=chain)

            tax_ids = query.one().taxonomy_id
            if tax_ids is None:
                if simplify:
                    return None
                return []

            tax_ids = [int(tax_id) for tax_id in tax_ids.split(',')]

        # with self.session() as session:
        #     query = session.query(mod.SpeciesMapping.species_id).\
        #         filter(mod.SpeciesMapping.species_mapping_id.in_(tax_ids))
        #     species_ids = [result.species_id for result in query]
        with self.session() as session:
            query = session.query(mod.TaxidSpeciesDomain.species_taxid).\
            filter(mod.TaxidSpeciesDomain.taxonomy_id.in_(tax_ids))
            species_ids = [result.species_taxid for result in query]

            if simplify:
                if len(species_ids) > 1:
                    return SYNTHEIC[0]
                if not species_ids:
                    return None
                return species_ids[0]

            if len(species_ids) != len(tax_ids):
                missing = set(species_ids).symmetric_difference(set(tax_ids))
                missing = [str(tax_id) for tax_id in missing]
                raise UnknownTaxonomyException("Missing tax ids %s",
                                               ','.join(missing))

            return sorted(species_ids)


class BasePairQueries(Base):
    """This is a class to deal with getting information about basepairs from
    the database. We store some useful but complex queries in here as methods
    for ease of use.
    """

    def representative(self, pdb, chain, model=1, count=False,
                       range_cutoff=None, near=False, family=None,
                       sym_op='1_555'):
        """This gets all forward interactions within a chain. This means we get
        only the interactions which are either non-symmetric (tWH) or the
        symmetric ones where unit 1 id < unit 2 id.

        :pdb: The pdb to search.
        :chain: The chain to search in.
        :count: If we should return the count or the interactions.
        :near: Should we count nears.
        :family: The family to limit to.
        :returns: A count or the interactions themself.
        """
        inter = mod.UnitPairsInteractions

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chain, near=near,
                                          family=family, model=model,
                                          sym_op=sym_op)

            query = query.filter(u1.chain == u2.chain)

            if range_cutoff is not None:
                query = query.filter(inter.f_crossing >= range_cutoff)

            if count:
                return query.count()
            return [result for result in query]

    def cross_chain(self, pdb, chain, other_chain=None, count=False,
                    near=False, family=None, model=1, sym_op='1_555'):
        """This method gets all interactions between on chain and another
        chain. If no other chain is specified then we go to any other chain.
        This will get all interactions

        :pdb: The pdb to search.
        :chain: The first chain.
        :other_chain: The second chain(s) if any.
        :count: If we should return a count or all interactions
        :near: True if we should count nears or not
        :family: The family to limit to.
        :returns: The count or a list of bps.
        """

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chain, near=near,
                                          family=family, symmetry=False,
                                          model=model, sym_op=sym_op)
            query = query.filter(u2.chain != u1.chain)

            if other_chain is not None:
                if isinstance(other_chain, list):
                    query = query.filter(u2.chain.in_(other_chain))
                else:
                    query = query.filter(u2.chain == other_chain)

            if count:
                return query.count()
            return [result for result in query]

    def between(self, pdb, chains, near=False, count=False, family=None,
                model=1, sym_op='1_555'):
        """This is a method to get the interactions between a list of chains.
        This is like a cross chain interaction but more general. Cross chain
        can be from A to B or C but not between B and C. This allows for
        interactions between A, B and C.

        :pdb: The pdb id.
        :chains: A list of chains.
        :near: A boolean to control showing nears.
        :count: A boolean to return a count or the list.
        :returns: A count or the list of interactions.
        """

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chains, near=near,
                                          family=family, model=model,
                                          sym_op=sym_op)
            query = query.\
                filter(u2.chain != u1.chain).\
                filter(u2.chain.in_(chains))

            if count:
                return query.count()
            return [result for result in query]

    def __base__(self, session, pdb, chain, symmetry=True, near=False,
                 family=None, model=1, sym_op='1_555'):
        """A method to build the base queries for this class.

        :session: A database session.
        :pdb: The pdb id.
        :chain: The chain(s) to query.
        :symmetry: If we should try to deduplicate symmetric basepairs.
        :near: If we should allow near interactions.
        :family: The family(ies) to limit the query to.
        :returns: A query.
        """

        u1 = aliased(mod.UnitInfo)
        u2 = aliased(mod.UnitInfo)
        bp = mod.BpFamilyInfo
        inter = mod.UnitPairsInteractions

        query = session.query(inter.unit_id_1, inter.unit_id_2, inter.f_lwbp).\
            join(u1, u1.unit_id == inter.unit_id_1).\
            join(u2, u2.unit_id == inter.unit_id_2).\
            join(bp, bp.bp_family_id == inter.f_lwbp).\
            filter(bp.is_forward == True).\
            filter(u1.model == u2.model).\
            filter(inter.pdb_id == pdb).\
            filter(u1.sym_op == u2.sym_op).\
            filter(u1.model == u2.model).\
            filter(u1.model == model).\
            filter(u1.sym_op == sym_op)

        if symmetry:
            query = query.filter((bp.is_symmetric == False) |
                                 ((bp.is_symmetric == True) & (u1.unit_id < u2.unit_id)))

        if not near:
            query = query.filter(bp.is_near == False)

        if family is not None:
            if isinstance(family, list):
                query = query.filter(inter.f_lwbp.in_(family))
            else:
                query = query.filter(inter.f_lwbp == family)

        if isinstance(chain, list):
            query = query.filter(u1.chain.in_(chain))
        elif chain:
            query = query.filter(u1.chain == chain)

        return u1, u2, query


class Correspondence(Base):
    def reference(self, pdb):
        """Get all correlated reference structures.
        """
        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo).filter_by(pdb2=pdb)
            return [ut.row2dict(result) for result in query]

    def mapping(self, ref, pdb):
        """Get the mapping from nucleotides in the reference to the nucleotides
        in the pdb.
        """
        mapping = {}
        with self.session() as session:
            query = session.query(mod.CorrespondenceNts).\
                join(mod.CorrespondenceInfo,
                     mod.CorrespondenceInfo.id ==
                     mod.CorrespondenceNts.correspondence_id).\
                filter(mod.CorrespondenceInfo.pdb1 == ref).\
                filter(mod.CorrespondenceInfo.pdb2 == pdb)
            for result in query:
                mapping[result.unit_id_1] = result.unit_id_2
                mapping[result.unit_id_2] = result.unit_id_1

        return mapping

    # def experimental_sequence_mapping(self, pdb, chain):
    #     with self.session() as session:
    #         query = session.query(ObsMap).\
    #             join(Obs, Obs.id == ObsMap.sequence_unit_id).\
    #             filter(ObsMap.pdb == pdb)
    #         return [ut.row2dict(result) for result in query]
