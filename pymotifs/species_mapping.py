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
    """A parser to take the response from a webrequest to the EBI and pull out
    the species assignment from the resulting XML document.
    """

    def get_species(self, taxon):
        """Given the taxon entry from a request extract the species id and
        scientific name. If there no assigned species or the taxon is above the
        species level, then we return (None, None).

        Parameters
        ----------
        taxon : ET.
            The taxon entry

        Returns
        -------
        assignment : (int, str)
            The assigned species id and scientific name for the species.
        """

        if taxon.attrib.get('rank') == 'species':
            return int(taxon.attrib['taxId']), taxon.attrib['scientificName']

        lineage = taxon.find('lineage')
        for level in lineage.findall('taxon'):
            if level.attrib.get('rank') == 'species':
                return (int(level.attrib['taxId']),
                        level.attrib['scientificName'])

        return None, None

    def parse(self, text):
        """Actually do the parsing. This method is separated out from the
        __call__ method for testing purposes.

        Parameters
        ----------
        text : str
            The text to parse

        Returns
        -------
        species : dict
            A dict with 'species_mapping_id', 'species_id' and 'species_name'
            fields. These fields contain the taxon id that was queried with,
            the species id and the species name respectively.
        """

        root = ET.fromstring(str(text))
        taxon = root.find('taxon')
        species_id, species_name = self.get_species(taxon)
        return {
            'species_mapping_id': int(taxon.attrib['taxId']),
            'species_id': species_id,
            'species_name': species_name
        }

    def __call__(self, response):
        """Parse the response text. The parsing is done by `Parser.parse`.

        Parameters
        ----------
        response : requests.Response
            The response to parse

        Returns
        -------
        species : dict
            A dict with 'species_mapping_id', 'species_id' and 'species_name'
            fields. These fields contain the taxon id that was queried with,
            the species id and the species name respectively.
        """
        return self.parse(response.text)


class Loader(core.SimpleLoader):
    """The loader to actually fetch and store the species assignments for all
    taxonomy ids.
    """

    """We allow this to merge data since sometimes we want to replace"""
    merge_data = True

    """The dependencies we require."""
    dependencies = set([ChainLoader])

    """URL pattern to use for taxon id data."""
    url = 'https://www.ebi.ac.uk/ena/data/view/Taxon:%s&display=xml'

    """The species id will not fit in the pdb id column used to mark"""
    mark = False

    def to_process(self, pdbs, **kwargs):
        """Determine the data to process. We do not process PDB ids in this
        loader, instead we get all taxon ids that are known about in the
        database. We want to look up the taxonomy for all of these. We do not
        work on the basis of on structure as there are both duplicates and one
        structure corresponds to several taxonomies as it may have several
        species.

        Parameters
        ----------
        pdbs : list
            A list of PDB ids to transform.

        Returns
        -------
        taxon_ids : list
            A list of taxon ids to lookup the taxonomy for.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.taxonomy_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                filter(mod.ChainInfo.taxonomy_id != None).\
                distinct().\
                order_by(mod.ChainInfo.taxonomy_id)

            ids = []
            for result in query:
                try:
                    ids.append(int(result.taxonomy_id))
                except:
                    self.logger.info("Bad taxonomy id %s", result.taxonomy_id)
            return ids

    def query(self, session, taxonomy_id):
        """Create a query to find the entires in species_mappings for the given
        taxonomy id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use
        taxonomy_id : int
            The taxonomy id to lookup

        Returns
        -------
        query : Query
            The query.
        """
        return session.query(mod.SpeciesMapping).\
            filter_by(species_mapping_id=taxonomy_id)

    def data(self, taxonomy_id, **kwargs):
        """Fetch the taxonomic assignment information for the given
        taxonomy_id. This will make a web request (using retries) to get the
        taxonomy information. It will then extract the species assignment and
        produce data for it.

        Parameters
        ----------
        taxonomy_id : int
            The taxonomy_id id to use.

        Returns
        -------
        species : pymotifs.models.SpeciesMapping
            The species mapping object to save.
        """

        helper = WebRequestHelper(parser=Parser())
        result = helper(self.url % taxonomy_id)
        if 'species_id' not in result:
            self.logger.warning("No species found for %i", taxonomy_id)
        return mod.SpeciesMapping(**result)
