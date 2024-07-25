"""
Map each taxid in chain_info to its species, domain, and lineage

This queries NCBI's taxononmy services to determine the species id for all
taxon ids assigned to chains.
Some chain organisms are species level and some are not.
This will determine the assigned species.
A few chains have assignments to the genus level, which this complains about.
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

    def get_rank_data(self, taxon, rank='species'):
        """
        Given the taxon entry from a request extract the tax id and
        scientific name.
        If there no matching rank or the taxon is above the
        species level, then we return (None, None).

        Parameters
        ----------
        taxon : ET.
            The taxon entry

        rank : species, superkingdom

        Returns
        -------
        assignment : (int, str)
            The assigned species id and scientific name for the species.
        """

        if taxon.attrib.get('rank') == rank:
            return (taxon.attrib['taxId'], taxon.attrib['scientificName'])

        lineage = taxon.find('lineage')
        for level in lineage.findall('taxon'):
            if level.attrib.get('rank') == rank:
                return (level.attrib['taxId'], level.attrib['scientificName'])

        return None, None


    def parse(self, text):
        """
        Actually do the parsing. This method is separated out from the
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

        dd = {}

        for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
            species_taxid, rank_value = self.get_rank_data(taxon,rank)

            if rank == 'superkingdom':
                if rank_value == 'Eukaryota':
                    rank_value = 'Eukarya'
                dd['domain'] = rank_value
            else:
                dd[rank] = rank_value

        dd['species_taxid'] = species_taxid     # taxid of the main species

        return dd


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

    """The species id will not fit in the pdb id column used to mark"""
    mark = False

    def to_process(self, pdbs, **kwargs):
        """
        Determine the taxonomy ids to process.
        Retrieve all taxonomy ids from the chain_info table.
        Only look at current pdb ids.

        Parameters
        ----------
        pdbs : list
            A list of PDB ids to focus on.

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

            needed_ids = set()
            for result in query:
                ids = result.taxonomy_id.split(",")
                for id in ids:
                    needed_ids.add(id)

        with self.session() as session:
            query = session.query(mod.TaxidSpeciesDomain.taxonomy_id).\
                distinct()

            known_ids = set()
            for result in query:
                known_ids.add(result.taxonomy_id)

        look_up = sorted(needed_ids - known_ids, key = lambda x : int(x))

        if len(look_up) == 0:
            raise core.Skip("No taxonomy ids to process")

        return look_up

    def query(self, session, taxid):
        """
        Create a query to find the entries in taxid_species_domain for the given
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

        return session.query(mod.TaxidSpeciesDomain).\
            filter_by(taxonomy_id=taxid)

    def data(self, taxid, **kwargs):
        """Fetch the taxonomic assignment information for the given
        taxid. This will make a web request (using retries) to get the
        taxonomy information. It will then extract the species assignment and
        produce data for it.

        Parameters
        ----------
        taxid : int
            The taxid id to use.

        Returns
        -------
        species : pymotifs.models.SpeciesMapping
            The species mapping object to save.
        """

        # manual mapping for a few tax ids that are not found at www.ebi.ac.uk
        # new ones are written to /usr/local/pipeline/hub-core/taxonomic_work_to_do.txt
        # find the ParentTaxId in the xml file and use that instead
        # TODO: look them up automatically at nih.gov; need a hole in the firewall for that
        manual = {}
        manual["11176"] = "2560319" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=11176&retmode=xml
        manual["11269"] = "3052505" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=11269&retmode=xml
        manual["11628"] = "3052317" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=11628&retmode=xml
        manual["121791"] = "3052225" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=121791&retmode=xml
        manual["1217710"] = "70346" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1217710&retmode=xml
        manual["1513237"] = "3052148" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1513237&retmode=xml
        manual["1933179"] = "3052381" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1933179&retmode=xml
        manual["1979162"] = "2560580" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1979162&retmode=xml
        manual["1980471"] = "3052480" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1980471&retmode=xml
        manual["2170080"] = "3052619" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2170080&retmode=xml
        manual["2170082"] = "3052625" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2170082&retmode=xml
        manual["2560547"] = "3052409" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2560547&retmode=xml
        manual["2560550"] = "3052410" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2560550&retmode=xml
        manual["2560743"] = "3052437" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2560743&retmode=xml
        manual["2734478"] = "3052686" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2734478&retmode=xml
        manual["1933298"] = "2560196" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1933298&retmode=xml

        if taxid in manual:
            new_taxid = manual[taxid]
        else:
            new_taxid = taxid

        """URL pattern to use for taxon id data."""
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s&retmode=xml'  # not working!
        url = 'https://www.ebi.ac.uk/ena/data/view/Taxon:%s&display=xml'

        try:
            helper = WebRequestHelper(parser=Parser())
            result = helper(url % new_taxid)
        except Exception as e:
            self.logger.warning("No taxonomic information found for %s" % taxid)
            manual_work = 'manual["%s"] = "" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s&retmode=xml' % (taxid,taxid)
            self.logger.info(manual_work)
            with open('taxonomic_work_to_do.txt','a') as f:
                f.write(manual_work + '\n')
            raise core.Skip("Unable to load some taxonomic information")

        result['taxonomy_id'] = taxid

        if not 'species_taxid' in result:
            self.logger.warning("No taxonomic information found for %s" % taxid)

            manual_work = 'manual["%s"] = "" # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s&retmode=xml' % (taxid,taxid)
            self.logger.info(manual_work)
            with open('taxonomic_work_to_do.txt','a') as f:
                f.write(manual_work + '\n')

            raise core.Skip("Unable to load some taxonomic information")
        # self.logger.info('taxid_species_domain stage works successfully')
        return mod.TaxidSpeciesDomain(**result)
