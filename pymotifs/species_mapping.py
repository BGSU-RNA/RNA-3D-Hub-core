"""Map from taxon id to species taxon id.

This queries NCBI's taxononmy services to determine the species id for all
taxon ids assigned to chains. Some chain organisms are species level and some
are not. This will determine the assigned species. A few chains have
assignments to the genus level, which this complains about.
"""

try:
    import lxml
except ImportError:
    pass

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import WebRequestHelper
from pymotifs.chains.info import Loader as ChainLoader


class Parser(object):
    def __call__(self, response):
        pass


class Loader(core.SimpleLoader):

    merge_data = True
    dependencies = set([ChainLoader])
    url = 'https://www.ebi.ac.uk/ena/data/view/Taxon:%s&display=xml'

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.ChainInfo.taxonomy_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                distinct()
            return [result.taxonomy_id for result in query]

    def query(self, session, taxonomy_id):
        return session.query(mod.SpeciesMapping).\
            filter_by(mod.SpeciesMapping.species_mapping_id == taxonomy_id)

    def data(self, taxonomy_id, **kwargs):
        helper = WebRequestHelper(parser=Parser())
        result = helper(self.url % taxonomy_id)
        return [mod.SpeciesMapping(**d) for d in result]
