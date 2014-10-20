import json
import logging

import core

from models import CorrespondenceNts as Corr
from models import CorrespondenceInfo as Info

from utils import WebRequestHelper
from correspondence.utils import StructureUtil

logger = logging.getLogger(__name__)

URL = 'http://localhost:8080/api/services/correlations'


class Parser(object):
    def __init__(self, additional=None):
        self.additional = additional or {}

    def __call__(self, text):
        data = []
        for row in json.loads(text):
            entry = {'unit1_id': row['unit1'], 'unit2_id': row['unit2']}
            entry.update(self.additional)
            data.append(entry)
        return data


class Loader(core.Loader):
    name = 'correspondence_nts'
    update_gap = False

    def __init__(self, config, maker):
        self.util = StructureUtil(maker)
        super(Loader, self).__init__(self, config, maker)

    def structure_data(self, chain, pdb):
        ids = self.util.unit_ids(pdb, chain)
        if not ids:
            logger.error("Could not load ids for PDB: %s, Chain: %s", pdb,
                         chain)
            logger.error("Was chain_breaks run prior to this")
            raise core.InvalidState("Could not find polymeric ids")

        sequence = self.util.polymer_sequences(pdb, chain)
        if not sequence:
            logger.error("Could not load sequence for PDB: %s, Chain: %s", pdb,
                         chain)
            logger.error("Was chain_breaks run prior to this")
            raise core.InvalidState("Could not find polymeric sequence")

        return {
            'ids': ids,
            'sequence': sequence,
            'pdb': pdb
        }

    def known(self, pdb):
        with self.session() as session:
            query = session.query(Info).filter(Info.pdb1 == pdb)
            return [(result.id, result.pdb2) for result in query]

    def correlate(self, correlation_id, ref, target):
        if not ref or not ref.get('ids') or not ref.get('sequence'):
            logger.error("Not given complete reference: %s", ref)
            raise core.InvalidState("Incomplete reference")

        if not target or not target.get('ids') or not target.get('sequence'):
            logger.error("Not given complete target: %s", target)
            raise core.InvalidState("Incomplete target")

        data = {
            'reference': ref['sequence'][0],
            'reference_ids': ref['ids'],
            'target': target['sequence'],
            'target_ids': target['ids']
        }

        helper = WebRequestHelper(method='post', parser=Parser())
        helper.parser.additional = {'correspondence_id': correlation_id}

        result = helper(URL, data=data, headers={'accept': 'application/json'})
        logger.info("Found %s correlations", len(result))
        return [Corr(**d) for d in result]

    def data(self, pdb, **kwargs):
        for corr_id, other in self.known(pdb):
            other_chain = self.util.longest_chain(other)
            logger.info("Using chain %s in %s", other_chain, other)

            pdb_chain = self.util.longest_chain(pdb)
            logger.info("Using chain %s in %s", pdb_chain, pdb)

            yield self.correlate(self.structure_data(pdb_chain, pdb),
                                 self.structure_data(other_chain, other))


if __name__ == '__main__':
    from utils import main
    main(Loader)
