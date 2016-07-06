from test import StageTest

from pymotifs.utils.naming import Namer


class ParentTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(ParentTest, self).setUp()
        self.group1 = {'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}
        self.group2 = {
            'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'D'}, {'id': 'E'}],
        }
        self.group3 = {'members': [{'id': 'D'}, {'id': 'E'}]}
        self.group4 = {'members': [{'id': 'E'}, {'id': 'F'}]}

    def test_it_can_compute_all_parents(self):
        val = self.loader.parents(self.group1,
                                  [self.group2, self.group3, self.group4])
        self.assertEquals(len(val), 1)

    def test_it_determines_the_size_of_each_overlap(self):
        val = self.loader.parents(self.group1,
                                  [self.group2, self.group3, self.group4])
        ans = [{'group': self.group2, 'intersection': set(['A', 'B'])}]
        self.assertEquals(ans, val)

    def test_it_can_find_several_parents(self):
        val = self.loader.parents(self.group2,
                                  [self.group1, self.group3, self.group4])
        self.assertEquals(len(val), 3)

    def test_it_gives_empty_list_for_no_parents(self):
        val = self.loader.parents(self.group1,
                                  [self.group3, self.group4])
        self.assertEquals([], val)


class OverlapTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(OverlapTest, self).setUp()
        self.group1 = {
            'members': [
                {'id': 'A'},
                {'id': 'B'},
                {'id': 'C'}
            ],
        }
        self.group2 = {
            'members': [
                {'id': 'A'},
                {'id': 'B'},
                {'id': 'D'},
                {'id': 'E'}
            ],
        }
        self.group3 = {
            'members': [
                {'id': 'D'},
                {'id': 'E'}
            ],
        }

    def test_it_can_determine_overlap(self):
        val = self.loader.overlap(self.group1, self.group2)
        ans = {'group': self.group2, 'intersection': set(['A', 'B'])}
        self.assertEqual(ans, val)

    def test_it_gives_empty_for_no_overlap(self):
        val = self.loader.overlap(self.group1, self.group3)
        ans = {}
        self.assertEqual(ans, val)

    def test_it_can_computes_the_min_overlap_size(self):
        overlap = self.loader.overlap(self.group1, self.group2)
        val = self.loader.overlap_size(self.group1, overlap)
        ans = 2.0 / 4.0
        self.assertEquals(ans, val)


class NewNamingTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(NewNamingTest, self).setUp()
        self.known = set(['10000', '20000'])

    def test_it_can_create_a_new_handler(self):
        naming = self.loader.new_name(1, self.known)
        self.assertEqual(5, len(naming['handle']))
        self.assertFalse(naming['handle'] in self.known)

    def test_sets_the_version_to_one(self):
        naming = self.loader.new_name(1, self.known)
        self.assertEqual(1, naming['version'])

    def test_sets_the_no_parent_comment(self):
        naming = self.loader.new_name(0, self.known)
        self.assertEqual('New id, no parents', naming['comment'])

    def test_sets_the_single_parent_comment(self):
        naming = self.loader.new_name(1, self.known)
        self.assertEqual('New id, 1 parent', naming['comment'])

    def test_sets_the_two_parent_comment(self):
        naming = self.loader.new_name(2, self.known)
        self.assertEqual('New id, 2 parents', naming['comment'])

    def test_sets_the_multiple_parent_comment(self):
        naming = self.loader.new_name(4, self.known)
        self.assertEqual('New id, > 2 parents', naming['comment'])


class SameNameTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(SameNameTest, self).setUp()
        self.group = {
            'handle': '00001',
            'version': 3,
            'comment': 'Something'
        }

    def test_it_copies_over_old_handle(self):
        name = self.loader.same_name(self.group)
        self.assertEqual('00001', name['handle'])

    def test_it_keeps_old_version(self):
        name = self.loader.same_name(self.group)
        self.assertEqual(3, name['version'])

    def test_it_adds_exact_comment(self):
        name = self.loader.same_name(self.group)
        self.assertEqual('Exact match', name['comment'])


class UpdatedNameTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(UpdatedNameTest, self).setUp()
        self.group = {
            'handle': '03001',
            'version': 3,
            'comment': 'Another comment'
        }

    def test_it_copies_over_old_handle(self):
        name = self.loader.updated_name(self.group, 1)
        self.assertEqual('03001', name['handle'])

    def test_it_bumps_old_version(self):
        name = self.loader.updated_name(self.group, 2)
        self.assertEqual(4, name['version'])

    def test_it_adds_comment_for_one_parent(self):
        name = self.loader.updated_name(self.group, 1)
        self.assertEqual('Updated, 1 parent', name['comment'])

    def test_it_adds_comment_for_two_parents(self):
        name = self.loader.updated_name(self.group, 2)
        self.assertEqual('Updated, 2 parents', name['comment'])


class OneParentTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(OneParentTest, self).setUp()
        self.group = {'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}
        self.exact = {
            'group': dict(self.group),
            'intersection': set(['A', 'B', 'C'])
        }
        self.exact['group']['name'] = {
            'handle': '10000',
            'comment': 'A',
            'version': 3
        }
        self.large = {
            'group': {
                'members': [{'id': 'A'}, {'id': 'C'}],
                'name': {'handle': '20000', 'comment': 'A', 'version': 3}
            },
            'intersection': set(['A', 'C'])
        }
        self.small = {
            'group': {
                'members': [{'id': 'D'}, {'id': 'C'}, {'id': 'E'}],
                'name': {'handle': '30000', 'comment': 'S', 'version': 10},
            },
            'intersection': set(['C'])
        }

    def test_it_uses_same_name_with_exact_overlap(self):
        val = self.loader.one_parent(self.group, self.exact, set([]))
        ans = {'handle': '10000', 'comment': 'Exact match', 'version': 3,
               'type': 'exact'}
        self.assertEqual(ans, val)

    def test_it_use_updated_name_with_large_overlap(self):
        val = self.loader.one_parent(self.group, self.large, set([]))
        ans = {'handle': '20000', 'comment': 'Updated, 1 parent', 'version': 4,
               'type': 'updated'}
        self.assertEqual(ans, val)

    def test_it_gives_new_name_if_overlap_to_old_is_small(self):
        name = self.loader.one_parent(self.group, self.small, set(['30000']))
        self.assertEqual(1, name['version'])
        self.assertEqual('New id, 1 parent', name['comment'])
        self.assertTrue(name['handle'] != '30000')

    def test_it_gves_new_name_if_overlap_to_new_is_small(self):
        group = dict(self.small['group'])
        overlap = dict(self.small)
        overlap['group']['members'] = dict(self.group)
        name = self.loader.one_parent(group, overlap, set(['30000']))
        self.assertEqual(1, name['version'])
        self.assertEqual('New id, 1 parent', name['comment'])
        self.assertTrue(name['handle'] != '30000')


class TwoParentTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(TwoParentTest, self).setUp()
        self.group = {'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}
        self.large = {
            'group': {
                'members': [{'id': 'A'}, {'id': 'C'}],
                'name': {'handle': '20000', 'comment': 'A', 'version': 3}
            },
            'intersection': set(['A', 'C'])
        }
        self.small = {
            'group': {
                'members': [{'id': 'D'}, {'id': 'C'}, {'id': 'E'}],
                'name': {'handle': '30000', 'comment': 'S', 'version': 10},
            },
            'intersection': set(['C'])
        }

    def test_it_updates_if_one_old_has_large_overlap(self):
        parents = [dict(self.large), dict(self.small)]
        known = set(['20000', '30000'])
        val = self.loader.two_parents(self.group, parents, known)
        ans = {
            'handle': '20000',
            'comment': 'Updated, 2 parents',
            'version': 4,
            'type': 'updated'
        }
        self.assertEqual(ans, val)

    def test_it_uses_new_if_both_have_small_overlap(self):
        parents = [dict(self.small), dict(self.small)]
        known = set(['20000', '30000'])
        name = self.loader.two_parents(self.group, parents, known)
        self.assertEquals(1, name['version'])
        self.assertTrue(name['handle'] not in known)
        self.assertEquals('New id, 2 parents', name['comment'])


class ManyParentTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(ManyParentTest, self).setUp()
        self.group = {'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}
        self.large = {
            'group': {
                'members': [{'id': 'A'}, {'id': 'C'}],
                'name': {'handle': '20000', 'comment': 'A', 'version': 3}
            },
            'intersection': set(['A', 'C'])
        }
        self.small = {
            'group': {
                'members': [{'id': 'D'}, {'id': 'C'}, {'id': 'E'}],
                'name': {'handle': '30000', 'comment': 'S', 'version': 10},
            },
            'intersection': set(['C'])
        }
        self.known = set(['20000', '30000'])

    def test_it_creates_new_name(self):
        parents = [dict(self.small), dict(self.small), dict(self.small)]
        name = self.loader.many_parents(self.group, parents, self.known)
        self.assertEquals(1, name['version'])
        self.assertTrue(name['handle'] not in self.known)
        self.assertEquals('New id, > 2 parents', name['comment'])

    def test_it_ignores_overlap_to_old(self):
        parents = [dict(self.large), dict(self.small), dict(self.small)]
        name = self.loader.many_parents(self.group, parents, self.known)
        self.assertEquals(1, name['version'])
        self.assertTrue(name['handle'] not in self.known)
        self.assertEquals('New id, > 2 parents', name['comment'])


class DataStructureTest(StageTest):
    loader_class = Namer

    def test_it_produces_the_correct_datastructure(self):
        groups = [{'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}]
        parents = [
            {
                'members': [{'id': 'A'}, {'id': 'C'}],
                'name': {'handle': '00001', 'version': 1}
            },
            {
                'members': [{'id': 'B'}, {'id': 'D'}, {'id': 'E'}],
                'name': {'handle': '20000', 'version': 3}
            }
        ]
        handles = set()
        val = self.loader(groups, parents, handles)
        ans = [{
            'comment': 'Updated, 2 parents',
            'parents': [
                {'name': {'version': 1, 'handle': '00001'},
                 'members': [{'id': 'A'}, {'id': 'C'}]},
                {'name': {'version': 3, 'handle': '20000'},
                 'members': [{'id': 'B'}, {'id': 'D'}, {'id': 'E'}]}],
            'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}],
            'name': {
                'handle': '00001',
                'version': 2,
                'type': 'updated',
            },
        }]
        self.assertEquals(ans, val)


# class RealisticTest(StageTest):
#     loader_class = Namer

#     def load_groups(self, release_id, include_names=True):
#         loader = ClassLoader(self.loader.config, self.loader.session.maker)
#         loaded = loader.load_release(release_id)
#         if include_names:
#             return loaded

#         for group in loaded:
#             group.pop('name')

#         return loaded

#     def setUp(self):
#         super(RealisticTest, self).setUp()
#         # Load the 0.16 release as it has many differences from 0.15
#         # http://rna.bgsu.edu/rna3dhub/nrlist/release/0.16/4.0A
#         self.target = self.load_groups('0.16', include_names=False)
#         self.parent = self.load_groups('0.15')
#         self.handles = set([])
#         self.named = self.loader(self.target, self.parent, self.handles)

#     def test_loaded_correct_number_of_groups(self):
#         self.assertEquals(334, len(self.named))

#     def test_it_can_get_all_new_names(self):
#         new_ids = filter(lambda g: g['type'] == 'new', self.named)
#         self.assertEquals(38, len(new_ids))

#     def test_it_can_update_all_groups(self):
#         updated = filter(lambda g: g['type'] == 'updated', self.named)
#         self.assertEquals(19, len(updated))

#     def test_it_correctly_updates(self):
#         updated = filter(lambda g: g['type'] == 'updated', self.named)
#         val = set(u['handle'] for u in updated)
#         ans = set(['12866', '33284', '35786', '43978', '18567',
#                    '70780', '30156', '80516', '40846', '32957',
#                    '46550', '91662', '82024', '17595', '89897',
#                    '39765', '37971', '41266', '25307'])
#         self.assertEqual(ans, val)

#     def test_it_can_find_all_exact_matches(self):
#         exact = filter(lambda g: g['type'] == 'exact', self.named)
#         self.assertEquals(277, len(exact))

#     def test_uses_correct_exact_handles(self):
#         val = sorted([g['handle'] for g in self.named if g['type'] == 'exact'])
#         ans = ['00642', '00834', '01065', '01208', '01504', '02048', '02301',
#                '02604', '02851', '02943', '03298', '03853', '04221', '05220',
#                '05908', '06681', '06778', '06871', '07085', '07210', '07374',
#                '08058', '08178', '08469', '09299', '09336', '09406', '09710',
#                '09744', '09919', '10365', '11098', '11265', '11340', '11464',
#                '11816', '11819', '12122', '12172', '12239', '12452', '13708',
#                '14020', '14197', '14414', '14526', '14913', '15089', '15733',
#                '15800', '16088', '16221', '16392', '16544', '16673', '16901',
#                '17130', '17198', '17384', '17511', '17574', '17591', '18344',
#                '18402', '18585', '18731', '19620', '20505', '20654', '20755',
#                '20846', '21101', '21412', '21817', '22966', '23158', '23257',
#                '23554', '23901', '23960', '24360', '24626', '24770', '25474',
#                '25990', '26024', '26026', '26065', '26351', '26679', '27351',
#                '27759', '28103', '29749', '30015', '30721', '31242', '31353',
#                '31517', '32011', '32214', '32244', '33239', '33913', '34633',
#                '35276', '35740', '36480', '36678', '36851', '38207', '38420',
#                '38700', '38702', '38798', '39310', '39983', '40299', '41576',
#                '42197', '42330', '42433', '42489', '42639', '42711', '42812',
#                '42991', '43323', '43810', '43865', '44042', '45325', '46204',
#                '47124', '47327', '47389', '47609', '48561', '48801', '48812',
#                '49495', '50215', '51873', '51877', '52372', '52511', '52733',
#                '52914', '53175', '53245', '53390', '53783', '53944', '54015',
#                '54441', '55098', '55157', '55160', '55278', '55462', '56140',
#                '56880', '56921', '57338', '58237', '58332', '58986', '59388',
#                '59546', '59580', '60283', '60542', '60620', '60745', '61569',
#                '61954', '62023', '62229', '63108', '63111', '63840', '64216',
#                '64636', '65398', '65411', '65520', '65527', '65782', '66249',
#                '66753', '66920', '67313', '67359', '67414', '67717', '68991',
#                '69460', '70388', '70618', '70729', '71431', '72606', '72617',
#                '72844', '73035', '73659', '73660', '73862', '74146', '74517',
#                '74969', '74973', '75586', '75939', '76076', '76083', '77315',
#                '77548', '78030', '78897', '78955', '79332', '79628', '79726',
#                '79810', '80232', '80517', '80518', '80544', '80869', '81563',
#                '81747', '81883', '82638', '83174', '84248', '84922', '85814',
#                '86147', '86343', '86551', '86592', '86811', '87307', '87602',
#                '87924', '88399', '88472', '88903', '89320', '89582', '89713',
#                '89878', '90154', '90169', '90868', '92609', '92621', '92723',
#                '92841', '92971', '93323', '93634', '93666', '94139', '95538',
#                '95970', '96517', '96825', '96866', '96943', '97609', '97850',
#                '98507', '98558', '98654', '99867']
#         self.assertEquals(ans, val)

#     def test_updates_handles(self):
#         # The handles should be modified to contain entries for each named
#         # group. It starts empty and then has 334 things added.
#         self.assertEqual(334, len(self.handles))
