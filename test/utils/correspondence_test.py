import string
import itertools as it

from test import StageTest

from pymotifs import models as mod
from pymotifs.utils.correspondence import Helper


class AlignedChainsTest(StageTest):
    loader_class = Helper

    def aligned(self, pdbs, **kwargs):
        return self.loader.aligned_chains(pdbs, use_names=True, **kwargs)

    def unknown_pdbs(self):
        comb = it.combinations(string.ascii_uppercase, 3)
        comb = it.imap(lambda c: ''.join(c), comb)
        return it.imap(lambda c: '0' + c, comb)

    def test_can_load_all_alignments_from_the_pdb(self):
        val = sorted(self.aligned(['1GID', '3T4B']).keys())
        assert ['1GID||A', '1GID||B', '3T4B||A'] == val

    def test_can_load_all_alignments_to_each(self):
        assert self.aligned(['1GID', '124D', '3T4B']) == {
            '3T4B||A': {'1GID||A': False, '1GID||B': False, '124D||B': False},
            '1GID||A': {'3T4B||A': False, '1GID||B': True, '124D||B': False},
            '1GID||B': {'3T4B||A': False, '1GID||A': True, '124D||B': False},
            '124D||B': {'1GID||A': False, '1GID||B': False, '124D||B': False}
        }

    def test_can_load_all_alignments_for_many_structures(self):
        pdbs = ['1GID', '3T4B'] * 4000
        pdbs.extend(self.unknown_pdbs())
        assert self.aligned(pdbs) == {
            '3T4B||A': {'1GID||A': False, '1GID||B': False},
            '1GID||A': {'3T4B||A': False, '1GID||B': True},
            '1GID||B': {'3T4B||A': False, '1GID||A': True},
        }

    def test_can_load_only_good_alignments(self):
        val = self.aligned(['1GID', '124D', '3T4B'], good=True)
        ans = {'1GID||A': {'1GID||B': True}, '1GID||B': {'1GID||A': True}}
        assert val == ans

    def test_can_load_only_bad_alignments(self):
        val = self.aligned(['1GID', '124D', '3T4B'], good=False)
        ans = {
            '3T4B||A': {'1GID||A': False, '1GID||B': False, '124D||B': False},
            '1GID||A': {'3T4B||A': False, '124D||B': False},
            '1GID||B': {'3T4B||A': False, '124D||B': False}
        }
        assert val == ans


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
        assert self.corr('1FCW', '4V42') == [
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

    def test_gives_nothing_if_no_corresponding(self):
        assert self.loader.chains('1DUH', '1KOG') == []


class CorrespondingPdbTest(StageTest):
    loader_class = Helper

    def test_can_get_pdbs_that_have_correspondence(self):
        assert self.loader.pdbs('1FCW') == ['1FCW', '4V42']


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
        assert self.ordering('4V9Q', 'BX', '4V9Q', 'DX') == {
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

    def test_gives_nothing_for_not_correspondence(self):
        assert self.ordering('1KOG', 'I', '1DUH', 'A') == {}


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
        assert self.mapping('4V9Q', 'BX', '4V9Q', 'DX') == {
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

    def test_it_gives_empty_if_no_mapping(self):
        assert self.mapping('1FJG', 'X', '1EKD', 'A') == {}
