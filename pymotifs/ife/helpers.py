import itertools as it
import collections as coll
import operator as op

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod
from pymotifs.constants import IFE_SEPERATOR
from pymotifs.constants import STRUCTURED_BP_COUNT
from pymotifs.utils import structures as st
from pymotifs.utils.sorting import total_ordering


class IfeLoader(core.Base):
    def load(self, pdb, chain):
        """This loads all information about a chain into a dictionary. This
        will load generic information about a chain, such as resolved, length,
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

        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id).\
                filter_by(pdb_id=data['pdb'], chain=data['chain'],
                          unit_type_id='rna')
            data['length'] = query.count()

        helper = st.BasePairQueries(self.session)
        rep = helper.representative
        data['internal'] = rep(data['pdb'], data['chain'], count=True,
                               family='cWW')

        return IfeChain(**data)

    def cross_chain_interactions(self, ifes):
        """Create a dictionary of the interactions between the listed chains.
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
        counter = helper.cross_chain
        for name1, name2 in pairs:
            count = counter(pdb, name1, name2, count=True, family='cWW')
            if name1 == name2:
                count = 0L
            interactions[name1][name2] = count

        return dict(interactions)

    def __call__(self, pdb):
        helper = st.Structure(self.session.maker)
        names = helper.rna_chains(pdb)
        ifes = [self.load(pdb, name) for name in names]
        return ifes, self.cross_chain_interactions(ifes)


@total_ordering
class IfeChain(object):
    def __init__(self, pdb=None, chain=None, internal=None, length=None,
                 full_length=None, db_id=None, bps=None):
        self.pdb = pdb
        self.chain = chain
        self.db_id = db_id
        self.internal = internal
        self.bps = bps
        self.length = length
        self.full_length = full_length

    @property
    def id(self):
        return '%s|%s' % (self.pdb, self.chain)

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
        if self.is_structured:
            chains = self.chains(structured=True)
            return IFE_SEPERATOR.join(c.id for c in chains)
        return self.chains()[0].id

    def chains(self, structured=None):
        fn = lambda c: True
        if structured is not None:
            fn = lambda c: c.is_structured == structured

        chains = it.ifilter(fn, self._chains)
        return sorted(chains, reverse=True)

    def add(self, chain):
        self._chains.add(chain)
        self.integral = max(self.integral, chain)
        if chain.is_structured:
            self.is_structured = True

    def merge(self, group):
        for chain in group.chains():
            self.add(chain)

    def __getattr__(self, key):
        if key in set(['pdb', 'internal', 'full_length', 'length',
                       'completeness']):
            if self.integral:
                return getattr(self.integral, key)
            return None
        super(IfeGroup, self).__getattr__(key)

    def __len__(self):
        return len(self._chains)

    def __ne__(self, other):
        return not other or self.id != other.id

    def __eq__(self, other):
        return other and self.id == other.id

    def __lt__(self, other):
        return other and self.id < other.id

    def __repr__(self):
        return '<IfeGroup %r (%r)>' % (self.id, self.integral)
