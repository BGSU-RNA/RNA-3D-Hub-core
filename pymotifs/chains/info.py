from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper
from pymotifs.models import ChainInfo


class Loader(core.MassLoader):
    merge_data = True

    names = {
        'structureId': 'pdb_id',
        'chainId': 'chain_id',
        'classification': 'classification',
        'macromoleculeType': 'macromolecule_type',
        'classification': 'classification',
        'macromoleculeType': 'macromolecule_type',
        'entityId': 'entity_id',
        'sequence': 'sequence',
        'chainLength': 'chain_length',
        'db_id': 'db_id',
        'db_name': 'db_name',
        'molecularWeight': 'molecular_weight',
        'secondaryStructure': "secondary_structure",
        'entityMacromoleculeType': 'entity_macromolecule_type',
        'Ki': 'ki',
        'Kd': 'kd',
        'compound': 'compound',
        'source': 'source',
        'taxonomyId': 'taxonomy_id',
        'biologicalProcess': 'biological_process',
        'cellularComponent': 'cellular_component'
    }

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.helper = CustomReportHelper(fields=self.names)

    def rename(self, report):
        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report.get(key)
        renamed['id'] = report.get('id')
        return renamed

    def get_ids(self, reports):
        pdbs = [report['structureId'] for report in reports]
        chains = [report['chainId'] for report in reports]

        mapping = {}
        known_pdbs = set()
        with self.session() as session:
            query = session.query(ChainInfo.id,
                                  ChainInfo.pdb_id,
                                  ChainInfo.chain_id).\
                filter(ChainInfo.pdb_id.in_(pdbs)).\
                filter(ChainInfo.chain_id.in_(chains))

            for result in query:
                mapping[(result.pdb_id, result.chain_id)] = result.id
                known_pdbs.add(result.pdb_id)

        for report in reports:
            db_id = mapping.get((report['structureId'], report['chainId']))
            if not db_id and report['structureId'] in known_pdbs:
                self.logger.error("It seems a new chain %s was added to %s",
                                  report['chainId'], report['structureId'])
            if db_id:
                report['id'] = db_id

        return reports

    def data(self, pdbs, **kwargs):
        reports = self.helper(pdbs)
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
