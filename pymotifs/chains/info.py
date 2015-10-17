from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper
from pymotifs.models import ChainInfo

from pymotifs.pdbs import Loader as PdbLoader


class Loader(core.MassLoader):
    merge_data = True
    dependencies = set([PdbLoader])

    names = {
        'structureId': 'pdb_id',
        'chainId': 'chain_name',
        'classification': 'classification',
        'macromoleculeType': 'macromolecule_type',
        'entityId': 'entity_name',
        'sequence': 'sequence',
        'chainLength': 'chain_length',
        'source': 'source',
        'taxonomyId': 'taxonomy_id',
        'entityMacromoleculeType': 'entity_macromolecule_type',
        'compound': 'compound'
    }

    def rename(self, report):
        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report.get(key)
        renamed['chain_id'] = report.get('chain_id')
        return renamed

    def get_ids(self, reports):
        pdbs = [report['structureId'] for report in reports]
        chains = [report['chainId'] for report in reports]

        mapping = {}
        known_pdbs = set()
        with self.session() as session:
            query = session.query(ChainInfo.chain_id,
                                  ChainInfo.pdb_id,
                                  ChainInfo.chain_name).\
                filter(ChainInfo.pdb_id.in_(pdbs)).\
                filter(ChainInfo.chain_name.in_(chains))

            for result in query:
                mapping[(result.pdb_id, result.chain_name)] = result.chain_id
                known_pdbs.add(result.pdb_id)

        for report in reports:
            db_id = mapping.get((report['structureId'], report['chainId']))
            if not db_id and report['structureId'] in known_pdbs:
                self.logger.error("It seems a new chain %s was added to %s",
                                  report['chainId'], report['structureId'])

            if db_id:
                report['chain_id'] = db_id

        return reports

    def data(self, pdbs, **kwargs):
        helper = CustomReportHelper(fields=self.names)
        reports = helper(pdbs)
        data = []
        for report in self.get_ids(reports):
            renamed = self.rename(report)
            data.append(ChainInfo(**renamed))

        known = set(pdbs)
        seen = set(entry.pdb_id for entry in data)

        if len(seen) != len(known):
            missing = known - seen
            self.logger.error("Could not get info on all pdbs: %s",
                              ','.join(missing))

        return data
