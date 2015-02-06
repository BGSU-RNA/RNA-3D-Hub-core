from nose import SkipTest

from test import StageTest

from pymotifs.species_mapping import Loader


class ProcessResultsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ProcessResultsTest, self).setUp()
        self.record = [
            {u'Lineage': 'other sequences; artificial sequences',
             u'Division': 'Synthetic',
             u'ParentTaxId': '81077',
             u'PubDate': '1993/04/27 01:00:00',
             u'LineageEx': [
                 {u'ScientificName': 'other sequences',
                  u'TaxId': '28384', u'Rank': 'no rank'},
                 {u'ScientificName': 'artificial sequences',
                  u'TaxId': '81077', u'Rank': 'no rank'}
             ],
             u'CreateDate': '1995/02/27 09:24:00',
             u'TaxId': '32630',
             u'Rank': 'species',
             u'ScientificName': 'synthetic construct'
             },
            {u'Lineage': 'Viruses; ssRNA viruses; ssRNA positive-strand viruses, no DNA stage; Flaviviridae; Hepacivirus; Hepatitis C virus; Hepatitis C virus genotype 4; Hepatitis C virus subtype 4a',
             u'Division': 'Viruses',
             u'ParentTaxId': '31653',
             u'PubDate': '2005/11/12 18:00:14',
             u'LineageEx': [
                 {u'ScientificName': 'Viruses', u'TaxId': '10239', u'Rank': 'superkingdom'},
                 {u'ScientificName': 'ssRNA viruses', u'TaxId': '439488', u'Rank': 'no rank'},
                 {u'ScientificName': 'ssRNA positive-strand viruses, no DNA stage', u'TaxId': '35278', u'Rank': 'no rank'},
                 {u'ScientificName': 'Flaviviridae', u'TaxId': '11050', u'Rank': 'family'},
                 {u'ScientificName': 'Hepacivirus', u'TaxId': '11102', u'Rank': 'genus'},
                 {u'ScientificName': 'Hepatitis C virus', u'TaxId': '11103', u'Rank': 'species'},
                 {u'ScientificName': 'Hepatitis C virus genotype 4', u'TaxId': '33745', u'Rank': 'no rank'},
                 {u'ScientificName': 'Hepatitis C virus subtype 4a', u'TaxId': '31653', u'Rank': 'no rank'}],
             u'TaxId': '356418',
             u'Rank': 'no rank',
             u'ScientificName': 'Hepatitis C virus ED43'
             }
        ]

    def test_can_process_a_record(self):
        val = self.loader.get_species(self.record[0])
        ans = {
            'id': 32630,
            'species_id': 32630,
            'species_name': 'synthetic construct'
        }
        self.assertEquals(ans, val)

    def test_can_determine_species(self):
        val = self.loader.get_species(self.record[1])
        ans = {
            'id': 356418,
            'species_id': 11103,
            'species_name': 'Hepatitis C virus'
        }
        self.assertEquals(ans, val)


class SpeciesMappingTest(StageTest):
    loader_class = Loader

    def test_downloads_a_context(self):
        context = self.loader.context([32630, 356418])
        val = self.loader.taxon_entries(context)
        val.sort(key=lambda v: v['id'])
        ans = [
            {'id': 32630, 'species_id': 32630, 'species_name': 'synthetic construct'},
            {'id': 356418, 'species_id': 11103, 'species_name': 'Hepatitis C virus'}
        ]
        self.assertEquals(ans, val)

    def test_can_find_all_missing_ids(self):
        raise SkipTest()
