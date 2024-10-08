import itertools as it
import functools as ft
import collections as coll
import operator as op

from sqlalchemy.sql.expression import func

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod
from pymotifs.constants import IFE_SEPARATOR
from pymotifs.constants import STRUCTURED_BP_COUNT
from pymotifs.utils import structures as st
from pymotifs.utils.sorting import total_ordering


class IfeLoader(core.Base):

    def best_model(self, pdb, sym_op):
        """
        Determine what model to use for ifes. We will use the model with the
        most basepairs. It tiebreaks on model number, lower is better.

        :pdb: The pdb id to use.
        :sym_op: The symmetry operator to use.
        :returns: The model number to use.
        """

        with self.session() as session:
            query = session.query(mod.UnitInfo.model).\
                filter_by(pdb_id=pdb).\
                distinct()
            models = [result.model for result in query]
            if not models:
                raise core.InvalidState("No models found for %s" % pdb)
            if len(models) == 1:
                return models[0]

        helper = st.BasePairQueries(self.session)
        count = ft.partial(helper.representative, pdb, None, count=True,
                           sym_op=sym_op)
        models = [(count(model=model), -1 * model) for model in models]
        return -1 * max(models)[1]

    def sym_op(self, pdb):
        """
        Pick a symmetry operator to work with. It doesn't really matter
        which one we use since they all have the same interactions by
        definition. So we just get all distinct and take the first.
        That is good enough.

        :param str pdb: The pdb id to use.
        :returns: The name of the symmetry operator.
        """

        with self.session() as session:
            return session.query(mod.UnitInfo.sym_op).\
                filter_by(pdb_id=pdb).\
                distinct().\
                limit(1).\
                one().\
                sym_op

    def load(self, pdb, chain, model=1, sym_op='1_555'):
        """
        This loads all information about a chain into a dictionary.
        This will load generic information about a chain, such as resolved, length,
        database id, the source and information about basepairing. The
        basepairing information is loaded from `bps`.

        :pdb: The pdb to search.
        :chain: The chain to search.
        :returns: A dictionary with
        """

        with self.session() as session:
            query = session.query(
                mod.ChainInfo.chain_id.label('db_id'),
                mod.ChainInfo.chain_name.label('chain'),
                mod.ChainInfo.pdb_id.label('pdb'),
                mod.ChainInfo.chain_length.label('full_length'),
            ).filter_by(chain_name=chain, pdb_id=pdb)

            data = ut.result2dict(query.one())
            data['model'] = model

        # with self.session() as session:
        #     query = session.query(mod.UnitInfo.sym_op,
        #                           func.count(1).label('count'),
        #                           ).\
        #         filter_by(pdb_id=data['pdb']).\
        #         filter_by(chain=data['chain']).\
        #         filter_by(model=model).\
        #         group_by(mod.UnitInfo.sym_op).\
        #         limit(1)
        #     result = query.first()
        #     data['length'] = 0
        #     if result:
        #         data['length'] = result.count
        # data['length'] will also calculate water units, this query works for rna structures
        # even if the rna has water units.
        # however, it is not work for dna structures for some unknown reason.
        # thus, ekko directly make it equal to the 'full_length'.
        # Dr. Zirbel believes we should used the full_length instead of making different lengths for hybrid chains
        data['length'] = data['full_length']

        helper = st.BasePairQueries(self.session)
        rep = helper.representative
        data['internal'] = rep(data['pdb'], data['chain'], count=True,
                               family='cWW', sym_op=sym_op)
        data['bps'] = rep(data['pdb'], data['chain'], count=True,
                          model=model, sym_op=sym_op)

        return IfeChain(**data)

    def cross_chain_interactions(self, ifes, sym_op='1_555'):
        """
        Create a dictionary of the interactions between the listed chains.
        This will get only the counts.

        :chains: A list of chain dictionaries.
        :returns: A dictionary of like { 'A': { 'B': 10 }, 'B': { 'A': 10 } }.
        """

        if not ifes:
            raise core.InvalidState("No ifes to get interactions between")

        pdb = ifes[0].pdb
        helper = st.BasePairQueries(self.session)
        interactions = coll.defaultdict(dict)
        pairs = it.product((ife.chain for ife in ifes), repeat=2)
        counter = ft.partial(helper.cross_chain, pdb, count=True, family='cWW',
                             sym_op=sym_op)
        for name1, name2 in pairs:
            count = counter(name1, name2)
            if name1 == name2:
                count = 0
            interactions[name1][name2] = count

        return dict(interactions)


    def __call__(self, pdb):
        helper = st.Structure(self.session.maker)
        names = helper.na_chains(pdb)
        sym_op = self.sym_op(pdb)
        model = self.best_model(pdb, sym_op)
        load = ft.partial(self.load, pdb, model=model, sym_op=sym_op)
        ifes = [load(name) for name in names]
        return ifes, self.cross_chain_interactions(ifes, sym_op=sym_op)


@total_ordering
class IfeChain(object):
    def __init__(self, pdb=None, chain=None, internal=None, length=None,
                 full_length=None, db_id=None, bps=None, model=None):
        self.pdb = pdb
        self.chain = chain
        self.db_id = db_id
        self.internal = internal
        self.bps = bps
        self.length = length
        self.full_length = full_length
        self.model = model

    @property
    def id(self):
        return '%s|%s|%s' % (self.pdb, self.model, self.chain)

    @property
    def is_structured(self):
        return self.internal >= STRUCTURED_BP_COUNT

    @property
    def completeness(self):
        if not self.full_length or not self.length:
            return None
        return self.length / self.full_length

    @property
    def bp_nt(self):
        return float(self.bps) / float(self.length)

    def __comp_attrs__(self):
        return op.attrgetter('bps', 'internal', 'length', 'full_length')

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        # this is how chains are sorted to decide which is listed first in the IFE
        if other is None:
            return False
        fn = self.__comp_attrs__()
        if fn(self) < fn(other):
            return True
        if fn(self) == fn(other):
            return self.chain > other.chain
        return False

    def __eq__(self, other):
        fn = self.__comp_attrs__()
        return other and fn(self) == fn(other) and self.chain == other.chain

    def __hash__(self):
        return hash(self.id)

    def __repr__(self):
        # this is the to-string method for IfeChain
        return '<IfeChain: %r (%r)>' % (self.id, self.internal)


@total_ordering
class IfeGroup(object):
    def __init__(self, *chains):
        self.is_structured = False
        self._chains = set()
        self.integral = None
        for chain in chains:
            self.add(chain)

    @property
    def id(self):
        # group id is where the ife id is stored
        # before 2024-06-20, the group in the "rest" case did not list all chains
        if self.is_structured:
            # originally a group is not structured
            # when you add a structured chain, it becomes structured
            # but then the list of chains is only the structured chains!
            # tRNA-mRNA is one use case; tRNA is structured, mRNA is not, don't make an IFE
            # however, 2OEU is a different case, A is structured, B is not, but make an IFE
            chains = self.chains(structured=True)
        else:
            chains = self.chains()

        # new code to cover the structured-structured case, the un-un case,
        # and some new structured-unstructured cases but not tRNA-mRNA
        all_chains = self.chains()
        if len(all_chains) == len(chains):
            # all chains are structured, or all are unstructured, or just one chain
            return IFE_SEPARATOR.join(c.id for c in chains)
        else:
            ife_id = IFE_SEPARATOR.join(c.id for c in all_chains)
            # pdb_id = ife_id.split('|')[0]
            # # cannot use self.logger.info here, using print instead
            # print("#### struct+unstruct new ife %s http://rna.bgsu.edu/rna3dhub/pdb/%s/2d" % (ife_id,pdb_id))

            # temporarily make code to remove old chains as being IFEs of their own
            # for c in all_chains:
            #     print("UPDATE ife_info SET new_style = 0 WHERE ife_id = '%s';" % (c.id))

            return ife_id

    def chains(self, structured=None):
        # return the chains ordered by the __lt__ method
        fn = lambda c: True
        if structured is not None:
            fn = lambda c: c.is_structured == structured

        chains = it.ifilter(fn, self._chains)
        return sorted(chains, reverse=True)

    def add(self, chain):
        self._chains.add(chain)

        # print("A self.integral is %s" % (self.integral))

        self.integral = max(self.integral, chain)

        # print("B self.integral is %s" % (self.integral))

        if chain.is_structured:
            self.is_structured = True

    def merge(self, group):
        for chain in group.chains():
            self.add(chain)

    def __getattr__(self, key):
        if key in set(['pdb', 'internal', 'full_length', 'length',
                       'completeness', 'bps', 'model']):
            if self.integral:
                return getattr(self.integral, key)
            return None
        raise AttributeError(key)

    def __len__(self):
        return len(self._chains)

    def __ne__(self, other):
        return not other or self.id != other.id

    def __eq__(self, other):
        return other and self.id == other.id

    def __lt__(self, other):
        return other and self.id < other.id

    def __repr__(self):
        # to-string method
        return '<IfeGroup: id %r (integral %r)>' % (self.id, self.integral)
