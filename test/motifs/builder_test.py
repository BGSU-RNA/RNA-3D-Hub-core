import pytest

import collections as coll

from test import StageTest

from pymotifs.motifs.builder import Known
from pymotifs.motifs.builder import Builder
from pymotifs.motifs.builder import Combiner


class CombinerTest(StageTest):
    loader_class = Combiner
    rel_id = '0.01'
    directory = 'test/files/motifs/IL_20120905_0000/'

    def empty(self):
        return coll.defaultdict(lambda: self.loader.empty_motif(self.rel_id))

    def data(self, method):
        fn = getattr(self.loader, method)
        return fn(self.directory, self.empty())

    def test_it_can_load_loops(self):
        data = self.data('loops')
        assert data['Group_268'] == {
            'release_id': self.rel_id,
            'members': [
                {'id': "IL_2GDI_002"},
                {'id': "IL_2GDI_006"},
                {'id': "IL_2HOJ_002"}
            ],
            'positions': [],
            'ordering': [],
            'signature': None,
        }

    def test_it_can_load_all_positions(self):
        data = self.data('positions')
        assert data['Group_268'] == {
            'release_id': self.rel_id,
            'members': [],
            'positions': [
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_20_U_', 'position': 1},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_21_G_', 'position': 2},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_38_C_', 'position': 3},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_39_U_', 'position': 4},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_40_G_', 'position': 5},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_41_A_', 'position': 6},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_42_G_', 'position': 7},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_43_A_', 'position': 8},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_44_A_', 'position': 9},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_45_A_', 'position': 10},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_20_U_', 'position': 1},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_21_G_', 'position': 2},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_38_C_', 'position': 3},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_39_U_', 'position': 4},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_40_G_', 'position': 5},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_41_A_', 'position': 6},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_42_G_', 'position': 7},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_43_A_', 'position': 8},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_44_A_', 'position': 9},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_45_A_', 'position': 10},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_20_U_', 'position': 1},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_21_G_', 'position': 2},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_38_C_', 'position': 3},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_39_U_', 'position': 4},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_40_G_', 'position': 5},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_41_A_', 'position': 6},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_42_G_', 'position': 7},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_43_A_', 'position': 8},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_44_A_', 'position': 9},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_45_A_', 'position': 10},
            ],
            'ordering': [],
            'signature': None,
        }

    def test_it_can_load_all_orderings(self):
        data = self.data('ordering')
        assert data['Group_268'] == {
            'release_id': self.rel_id,
            'members': [],
            'positions': [],
            'ordering': [
                {'loop_id': 'IL_2GDI_002', 'original_order': 1, 'similarity_order': 1},
                {'loop_id': 'IL_2GDI_006', 'original_order': 2, 'similarity_order': 2},
                {'loop_id': 'IL_2HOJ_002', 'original_order': 3, 'similarity_order': 3},
            ],
            'signature': None,
        }

    def test_it_can_load_all_signatures(self):
        data = self.data('signature')
        assert data['Group_268'] == {
            'release_id': self.rel_id,
            'members': [],
            'positions': [],
            'ordering': [],
            'signature': "10_cWW-cWW--tWH",
        }

    def test_it_can_load_and_merge_all_data(self):
        data = self.loader(self.rel_id, self.directory)
        assert data['Group_268'] == {
            'release_id': self.rel_id,
            'members': [
                {'id': "IL_2GDI_002"},
                {'id': "IL_2GDI_006"},
                {'id': "IL_2HOJ_002"}
            ],
            'positions': [
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_20_U_', 'position': 1},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_21_G_', 'position': 2},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_38_C_', 'position': 3},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_39_U_', 'position': 4},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_40_G_', 'position': 5},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_41_A_', 'position': 6},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_42_G_', 'position': 7},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_43_A_', 'position': 8},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_44_A_', 'position': 9},
                {'loop_id': 'IL_2GDI_002', 'unit_id': '2GDI_AU_1_X_45_A_', 'position': 10},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_20_U_', 'position': 1},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_21_G_', 'position': 2},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_38_C_', 'position': 3},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_39_U_', 'position': 4},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_40_G_', 'position': 5},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_41_A_', 'position': 6},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_42_G_', 'position': 7},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_43_A_', 'position': 8},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_44_A_', 'position': 9},
                {'loop_id': 'IL_2GDI_006', 'unit_id': '2GDI_AU_1_Y_45_A_', 'position': 10},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_20_U_', 'position': 1},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_21_G_', 'position': 2},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_38_C_', 'position': 3},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_39_U_', 'position': 4},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_40_G_', 'position': 5},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_41_A_', 'position': 6},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_42_G_', 'position': 7},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_43_A_', 'position': 8},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_44_A_', 'position': 9},
                {'loop_id': 'IL_2HOJ_002', 'unit_id': '2HOJ_AU_1_A_45_A_', 'position': 10},
            ],
            'ordering': [
                {'loop_id': 'IL_2GDI_002', 'original_order': 1, 'similarity_order': 1},
                {'loop_id': 'IL_2GDI_006', 'original_order': 2, 'similarity_order': 2},
                {'loop_id': 'IL_2HOJ_002', 'original_order': 3, 'similarity_order': 3},
            ],
            'signature': "10_cWW-cWW--tWH",
        }


class KnownTests(StageTest):
    loader_class = Known

    def test_it_can_load_all_loops(self):
        val = self.loader.loops('IL', '1.0')
        assert len(val) == 0

    def test_it_can_load_all_handles(self):
        val = self.loader.handles()
        assert len(val) == 0

    def test_it_can_load_all_names(self):
        val = self.loader.names('IL', '1.0')
        assert len(val) == 0
        # assert '' in val

    def test_it_can_load_all_motifs(self):
        val = self.loader.motifs('IL', '1.0')
        assert val == []


class BuilderTests(StageTest):
    loader_class = Builder
    directory = 'test/files/motifs/IL_20120905_0000/'

    def test_it_can_load_all_mutual_discrepancy(self):
        val = self.loader.mutual_discrepancy(self.directory)
        assert len(val) == 47992
        assert val[0] == {
            'loop_id_1': 'IL_1ZEV_001',
            'discrepancy': 0.4232,
            'loop_id_2': 'IL_4A1B_028',
        }

    @pytest.mark.skip()
    def test_it_can_compute_parent_counts(self):
        # val = self.loader.parent_counts(
        pass

    @pytest.mark.skip()
    def test_it_can_name_all_new_motifs(self):
        pass
        # data = self.loader('IL', '0.0', '0.1',

    @pytest.mark.skip()
    def test_it_assigns_the_release_ids(self):
        loops = []
        data = self.loader('IL', '0.1', '0.2', loops)
        assert data['release'] == '0.2'
        assert data['parent'] == '0.1'

    @pytest.mark.skip()
    def test_it_can_load_all_release_data(self):
        pass
