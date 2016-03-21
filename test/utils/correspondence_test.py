from test import StageTest

from pymotifs import models as mod
from pymotifs.utils.correspondence import Helper


class AlignedChainsTest(StageTest):
    loader_class = Helper

    def setUp(self):
        super(AlignedChainsTest, self).setUp()
        self.alignments = self.loader.aligned_chains(['1GID', '3T4B'])

    def test_can_load_all_alignments_from_the_pdb(self):
        self.assertEquals([60, 61, 172], sorted(self.alignments.keys()))

    def test_can_load_all_alignments_to_each(self):
        ans = {
            172L: {60L: False, 61L: False},
            60L: {172L: False, 61L: True},
            61L: {172L: False, 60L: True}
        }
        self.assertEquals(ans, self.alignments)

    def test_can_load_only_good_alignments(self):
        val = self.loader.aligned_chains(['1GID', '3T4B'], good=True)
        ans = {60L: {61L: True}, 61L: {60L: True}}
        self.assertEquals(ans, val)

    def test_can_load_only_bad_alignments(self):
        val = self.loader.aligned_chains(['1GID', '3T4B'], good=False)
        ans = {
            172L: {60L: False, 61L: False},
            60L: {172L: False},
            61L: {172L: False}
        }
        self.assertEquals(ans, val)


class ChainTest(StageTest):
    loader_class = Helper

    def corr(self, pdb1, pdb2):
        vals = self.loader.chains('1FCW', '4V42')
        corr = []
        for val in vals:
            entry = []
            for part in val[1:]:
                entry.append({'pdb': part['pdb'], 'name': part['name']})
            corr.append(tuple(entry))
        return corr

    def test_can_all_corresponding_chains(self):
        val = self.corr('1FCW', '4V42')
        ans = [
            ({'pdb': '1FCW', 'name': 'A'}, {'pdb': '4V42', 'name': 'AB'}),
            ({'pdb': '1FCW', 'name': 'B'}, {'pdb': '4V42', 'name': 'AB'}),
            ({'pdb': '1FCW', 'name': 'C'}, {'pdb': '4V42', 'name': 'AB'}),
            ({'pdb': '1FCW', 'name': 'D'}, {'pdb': '4V42', 'name': 'AB'}),
            ({'pdb': '1FCW', 'name': 'E'}, {'pdb': '4V42', 'name': 'AB'}),
            ({'pdb': '1FCW', 'name': 'A'}, {'pdb': '4V42', 'name': 'AC'}),
            ({'pdb': '1FCW', 'name': 'B'}, {'pdb': '4V42', 'name': 'AC'}),
            ({'pdb': '1FCW', 'name': 'C'}, {'pdb': '4V42', 'name': 'AC'}),
            ({'pdb': '1FCW', 'name': 'D'}, {'pdb': '4V42', 'name': 'AC'}),
            ({'pdb': '1FCW', 'name': 'E'}, {'pdb': '4V42', 'name': 'AC'}),
        ]
        self.assertEquals(val, ans)

    def test_gives_nothing_if_no_corresponding(self):
        self.assertEquals([], self.loader.chains('1DUH', '1KOG'))


class CorrespondingPdbTest(StageTest):
    loader_class = Helper

    def test_can_get_pdbs_that_have_correspondence(self):
        val = self.loader.pdbs('1FCW')
        self.assertEquals(['1FCW', '4V42'], val)


class OrderingTest(StageTest):
    loader_class = Helper

    def chain_id(self, pdb, name):
        with self.loader.session() as session:
            return session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=name).\
                one().\
                chain_id

    def corr_id(self, chain1, chain2):
        with self.loader.session() as session:
            return session.query(mod.CorrespondencePdbs).\
                filter_by(chain_id_1=chain1).\
                filter_by(chain_id_2=chain2).\
                first().\
                correspondence_id

    def ordering(self, pdb1, name1, pdb2, name2):
        c1 = self.chain_id(pdb1, name1)
        c2 = self.chain_id(pdb2, name2)
        corr_id = self.corr_id(c1, c2)
        return self.loader.ordering(corr_id, c1, c2)

    def test_can_get_ordering_alignments_between_chains(self):
        ans = {
            '4V9Q|1|BX|A|14': 0,
            '4V9Q|1|DX|A|14': 0,
            '4V9Q|1|BX|A|15': 1,
            '4V9Q|1|DX|A|15': 1,
            '4V9Q|1|BX|A|16': 2,
            '4V9Q|1|DX|A|16': 2,
            '4V9Q|1|BX|U|17': 3,
            '4V9Q|1|DX|U|17': 3,
            '4V9Q|1|BX|G|18': 4,
            '4V9Q|1|DX|G|18': 4,
        }
        self.assertEquals(ans, self.ordering('4V9Q', 'BX', '4V9Q', 'DX'))

    def test_gives_nothing_for_not_correspondence(self):
        self.assertEquals({}, self.ordering('1KOG', 'I', '1DUH', 'A'))


class MappingTest(StageTest):
    loader_class = Helper

    def corr_id(self, chain1, chain2):
        with self.loader.session() as session:
            return session.query(mod.CorrespondencePdbs).\
                filter_by(chain_id_1=chain1['id']).\
                filter_by(chain_id_2=chain2['id']).\
                first().\
                correspondence_id

    def chain(self, pdb, name):
        with self.loader.session() as session:
            chain_id = session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=name).\
                one().\
                chain_id
            return {'pdb': pdb, 'name': name, 'id': chain_id}

    def mapping(self, pdb1, name1, pdb2, name2):
        c1 = self.chain(pdb1, name1)
        c2 = self.chain(pdb2, name2)
        corr_id = self.corr_id(c1, c2)
        return self.loader.mapping(corr_id, c1, c2)

    def test_it_gives_mapping_between_chains(self):
        ans = {
            '4V9Q|1|BX|A|14': '4V9Q|1|DX|A|14',
            '4V9Q|1|BX|A|15': '4V9Q|1|DX|A|15',
            '4V9Q|1|BX|A|16': '4V9Q|1|DX|A|16',
            '4V9Q|1|BX|U|17': '4V9Q|1|DX|U|17',
            '4V9Q|1|BX|G|18': '4V9Q|1|DX|G|18',
            '4V9Q|1|DX|A|14': '4V9Q|1|BX|A|14',
            '4V9Q|1|DX|A|15': '4V9Q|1|BX|A|15',
            '4V9Q|1|DX|A|16': '4V9Q|1|BX|A|16',
            '4V9Q|1|DX|U|17': '4V9Q|1|BX|U|17',
            '4V9Q|1|DX|G|18': '4V9Q|1|BX|G|18',
        }
        self.assertEquals(ans, self.mapping('4V9Q', 'BX', '4V9Q', 'DX'))

    def test_it_gives_empty_if_no_mapping(self):
        self.assertEquals({}, self.mapping('1FJG', 'X', '1EKD', 'A'))
