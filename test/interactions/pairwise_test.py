from test import StageTest

from pymotifs.interactions.pairwise import Loader


class DetectingInteractionTypeTest(StageTest):
    loader_class = Loader

    def test_can_determine_if_bp(self):
        val = self.loader.interaction_type('cHS')
        ans = 'f_lwbp'
        self.assertEquals(val, ans)

    def test_can_determine_if_bph(self):
        val = self.loader.interaction_type('7BPh')
        ans = 'f_bphs'
        self.assertEquals(val, ans)

    def test_can_determine_if_base_ribose(self):
        val = self.loader.interaction_type('0BR')
        ans = 'f_brbs'
        self.assertEquals(val, ans)

    def test_can_determine_if_stack(self):
        val = self.loader.interaction_type('s55')
        ans = 'f_stacks'
        self.assertEquals(val, ans)

    def test_return_none_otherwise(self):
        val = self.loader.interaction_type('cH')
        self.assertTrue(val is None)

    def test_allows_water_as_bp(self):
        val = self.loader.interaction_type('wat')
        ans = 'f_lwbp'
        self.assertEquals(val, ans)

    def test_it_does_not_allow_perp(self):
        self.assertTrue(self.loader.interaction_type('perp') is None)
        self.assertTrue(self.loader.interaction_type('nperp') is None)

    def test_it_does_not_allow_rib(self):
        self.assertTrue(self.loader.interaction_type('rIB') is None)
        self.assertTrue(self.loader.interaction_type('nrIB') is None)


class ParsingACsvTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ParsingACsvTest, self).setUp()
        self.data = self.loader.\
            interactions_from_csv('files/interactions/1GID.csv', '1GID')

    def test_can_parse_all_enteris(self):
        val = len(self.data)
        ans = 1513
        self.assertEqual(ans, val)

    def test_it_converts_crossing_to_int(self):
        val = self.data[0].f_crossing
        ans = 0
        self.assertEqual(ans, val)

    def test_skips_invalid_entries(self):
        ints = [d.f_lwbp for d in self.data]
        self.assertFalse('perp' in ints)
