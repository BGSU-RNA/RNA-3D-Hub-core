"""Map from taxon id to species taxon id.

This queries NCBI's taxononmy services to determine the species id for all
taxon ids assigned to chains. Some chain organisms are species level and some
are not. This will determine the assigned species. A few chains have
assignments to the genus level, which this complains about.
"""

import xml.etree.ElementTree as ET

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import WebRequestHelper
from pymotifs.chains.info import Loader as ChainLoader


class Parser(object):
    def get_species(self, taxon):
        if taxon.attrib.get('rank') == 'species':
            return int(taxon.attrib['taxId']), taxon.attrib['scientificName']

        lineage = taxon.find('lineage')
        for level in lineage.findall('taxon'):
            if level.attrib.get('rank') == 'species':
                return int(level.attrib['taxId']), level.attrib['scientificName']

        return None, None

    def parse(self, text):
        root = ET.fromstring(text)
        taxon = root.find('taxon')
        species_id, species_name = self.get_species(taxon)
        return {
            'species_mapping_id': int(taxon.attrib['taxId']),
            'species_id': species_id,
            'species_name': species_name
        }

    def __call__(self, response):
        return self.parse(response.text)


class Loader(core.SimpleLoader):

    merge_data = True
    dependencies = set([ChainLoader])
    url = 'https://www.ebi.ac.uk/ena/data/view/Taxon:%s&display=xml'
    mark = False

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.ChainInfo.taxonomy_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                filter(mod.ChainInfo.taxonomy_id != None).\
                distinct().\
                order_by(mod.ChainInfo.taxonomy_id)
            return [int(result.taxonomy_id) for result in query]

    def query(self, session, taxonomy_id):
        return session.query(mod.SpeciesMapping).\
            filter_by(species_mapping_id=taxonomy_id)

    def data(self, taxonomy_id, **kwargs):
        helper = WebRequestHelper(parser=Parser())
        result = helper(self.url % taxonomy_id)
        if result.species_id is None:
            self.logger.warning("No species found for %i", taxonomy_id)
        return mod.SpeciesMapping(**result)
