from test import StageTest

from pymotifs.correspondence.summary import Loader


class CheckingSmallAlignmentsTest(StageTest):
    loader_class = Loader

    def test_only_compares_alignments_of_same_size(self):
        info = {'aligned_count': 18, 'mismatch_count': 0}
        self.assertFalse(self.loader.good_alignment(info, 17, 18))

    def test_requires_identical_sequences(self):
        info = {'aligned_count': 18, 'mismatch_count': 0}
        self.assertTrue(self.loader.good_alignment(info, 18, 18))

    def test_fails_if_any_mismatches(self):
        info = {'aligned_count': 18, 'mismatch_count': 1}
        self.assertFalse(self.loader.good_alignment(info, 18, 18))

    def test_fails_if_none_aligned(self):
        info = {'aligned_count': 0, 'mismatch_count': 0}
        self.assertFalse(self.loader.good_alignment(info, 18, 18))

    def test_fails_if_all_aligned_but_differ_from_size(self):
        info = {'aligned_count': 18, 'mismatch_count': 1}
        self.assertFalse(self.loader.good_alignment(info, 17, 18))


class CheckingMediumAlignmentsTest(StageTest):
    loader_class = Loader

    def test_fails_with_too_many_mistmatches(self):
        info = {'aligned_count': 79, 'match_count': 70, 'mismatch_count': 9}
        self.assertFalse(self.loader.good_alignment(info, 79, 79))

    def test_passes_with_mismatch_less_than_5(self):
        info = {'aligned_count': 79, 'match_count': 75, 'mismatch_count': 4}
        self.assertTrue(self.loader.good_alignment(info, 79, 79))

    def test_fails_if_none_aligned(self):
        info = {'aligned_count': 0, 'match_count': 0, 'mismatch_count': 0}
        self.assertFalse(self.loader.good_alignment(info, 78, 78))

    def test_fails_with_realistic_example(self):
        info = {'aligned_count': 25, 'match_count': 23, 'mismatch_count': 5}
        self.assertFalse(self.loader.good_alignment(info, 28, 26))


class CheckingLargeAlignmentsTest(StageTest):
    loader_class = Loader

    def test_fails_if_none_aligned(self):
        info = {'aligned_count': 0, 'mismatch_count': 0}
        self.assertFalse(self.loader.good_alignment(info, 100, 100))

    def test_passes_if_no_mismatches(self):
        info = {'aligned_count': 100, 'mismatch_count': 0}
        self.assertTrue(self.loader.good_alignment(info, 100, 100))

    def test_fails_with_less_than_95_percent_similarity(self):
        info = {'aligned_count': 100, 'match_count': 94, 'mismatch_count': 6}
        self.assertFalse(self.loader.good_alignment(info, 100, 100))

    def test_fails_with_less_than_95_percent_similarity_ignoring_longest(self):
        info = {'aligned_count': 100, 'match_count': 94, 'mismatch_count': 6}
        self.assertFalse(self.loader.good_alignment(info, 100, 200))

    def test_passes_with_95_percent_similarity_ignoring_longest(self):
        info = {'aligned_count': 100, 'match_count': 95, 'mismatch_count': 5}
        self.assertTrue(self.loader.good_alignment(info, 100, 200))

    def test_it_passes_realistic_example(self):
        info = {'length': 1522, 'aligned_count': 1511, 'match_count': 1510,
                'mismatch_count': 1}
        self.assertTrue(self.loader.good_alignment(info, 1511, 1522))
