from Bio import Entrez

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.MassLoader):
    """This is a loader to get the species information for all structures. What
    this is extract all taxon ids for all chains that we do not yet have a
    species level information for then use NCIB's Entrez EPost to store all
    ids. We then later download all the information in large batches using the
    history.  We must do this because Entrez has a policy on the number and
    speed of requests. If we send too much too fast we run the risk of being
    banned. The current setup for this stage is unlikely to cause such problems
    as we can download all mapping information in a single request. In
    addition, we use BioPython's Entrez package which is supposed to limit
    requests to requested speed.
    """

    """The maximum number of ids to post to one WebContext."""
    context_limit = 3000

    """The number of ids to download in one request at once."""
    download_limit = 200

    merge_data = True

    allow_no_data = True

    dependencies = set([ChainLoader])

    def structure_taxon_ids(self, *args, **kwargs):
        """This transform is a bit special. We actually ignore the input to the
        function and simply find all chain entries that we don't know the
        species for. We then make a request to NCBI using history to store all
        the taxon ids. Later on we retrieve all ids using the history. We
        upload 3000 ids per WebContext and then download them in chunks later.

        :returns: A list of WebContext ids.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.taxonomy_id).\
                outerjoin(mod.SpeciesMapping,
                          mod.SpeciesMapping.species_mapping_id == mod.ChainInfo.taxonomy_id).\
                filter(mod.SpeciesMapping.species_mapping_id == None)

            ids = []
            for result in query:
                if result[0]:
                    for part in result[0].split(','):
                        ids.extend(part.split('#'))
            return set(ids)

    def contexts(self, ids):
        Entrez.email = self.config['email']
        grouped = utils.grouper(self.context_limit, ids)
        return {
            'ids': ids,
            'contexts': [self.context(chunk) for chunk in grouped]
        }

    def context(self, ids):
        """Given a list of ids we post them to entrez and return a dictonary
        with the important information such as web context (stored as 'id'),
        the query key (stored as 'key') and the number of requested ids (stored
        as 'size').

        :ids: A list of ids to post.
        :returns: A dictonary with 'id', 'key' and 'size' keys.
        """

        Entrez.email = self.config['email']
        str_ids = [str(id) for id in ids]
        handle = Entrez.epost("taxonomy", id=','.join(str_ids))
        data = Entrez.read(handle)
        return {
            'id': data['WebEnv'],
            'key': data['QueryKey'],
            'size': len(ids),
            'ids': ids
        }

    def taxon_entries(self, context):
        """This takes a context information and downloads all requested taxon
        information. It will produce dictionaries stating the species level
        mapping for each

        :context: A context dictionary.
        :returns: A list of species dictionaries.
        """

        data = []
        for start in range(0, context['size'], self.download_limit):
            handle = Entrez.efetch(db="taxonomy",
                                   retstart=start,
                                   webenv=context['id'],
                                   query_key=context['key'],
                                   retmax=self.download_limit)
            records = Entrez.read(handle)
            for record in records:
                data.append(self.get_species(record))

        if len(data) != context['size']:
            known = [datum.get('species_mapping_id') for datum in data]
            missing = set(context['ids']) - set(filter(None, known))
            self.logger.warn("Did not get all ids for context: %s",
                             context['id'])
            self.logger.warn("Missing ids: %s", ','.join(missing))
        return data

    def get_species(self, record):
        """This process a single species record to produce a dictionary mapping
        the record id to the species id.

        :record: A dictionary produced by Entrez.read.
        :returns: A dictionary with the keys, 'id' the record taxon id and
        'species_id', representing the species id for this record.
        """

        tax_id = int(record['TaxId'])
        species = None

        if record['Rank'] == 'species':
            species = record
        else:
            species = filter(lambda e: e['Rank'] == 'species',
                             record['LineageEx'])

        if not species:
            self.logger.warn("Could not determine species for %s", tax_id)
            return {}

        if isinstance(species, list):
            species = species[0]

        return {
            'species_mapping_id': tax_id,
            'species_id': int(species['TaxId']),
            'species_name': species['ScientificName']
        }

    def data(self, *args, **kwargs):
        """Get the SpeciesMapping for all contexts. If we can't get at least the
        expected number of objects warnings will be raised.

        :contexts: A dictionary produced by `transform`.
        :returns: A list of SpeciesMapping objects representing the species
        info.
        """

        ids = self.structure_taxon_ids()
        contexts = self.contexts(ids)
        data = []
        for context in contexts['contexts']:
            for entry in self.taxon_entries(context):
                if entry:
                    data.append(mod.SpeciesMapping(**entry))

        data = set(data)
        if len(data) < len(contexts['ids']):
            known = set([datum.species_mapping_id for datum in data])
            missing = set(context['ids']) - known
            self.logger.warn("Could not download all requested taxon ids")
            self.logger.warn("Missing: %s", ','.join(missing))

        return data
