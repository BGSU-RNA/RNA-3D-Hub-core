import json
import logging

from rnastructure.tertiary.cif import CIF

import core as core

from models import CorrespondenceNts as Corr
from models import CorrespondenceInfo as Info

import utils as utils
from utils.structures import NR as NrUtil
from utils.structures import Structure as StructureUtil

logger = logging.getLogger(__name__)

URL = 'http://localhost:8080/api/services/correlations'

PLACEHOLDER = '-'


class Parser(object):
    def __init__(self, additional=None):
        self.additional = additional or {}

    def __call__(self, response):
        data = []
        for row in json.loads(response.text):
            if row['unit1'] == PLACEHOLDER or row['unit2'] == PLACEHOLDER:
                continue
            entry = {'unit1_id': row['unit1'], 'unit2_id': row['unit2']}
            entry.update(self.additional)
            data.append(entry)
        return data


class Loader(core.Loader):
    name = 'correspondence_nts'
    update_gap = False

    def __init__(self, config, maker):
        self.nr_util = NrUtil(maker)
        self.st_util = StructureUtil(maker)
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def structure_data(self, chain, pdb):
        with open(self.finder(pdb), 'rb') as raw:
            structure = CIF(raw)

        ids = []
        sequence = []
        # TODO: Write out more generic unit ids
        chain_block = structure.chain('1_555', 1, chain)
        for (seq, _, unit_id) in chain_block.experimental_sequence_mapping():
            if unit_id is None:
                ids.append(PLACEHOLDER)
            else:
                ids.append(unit_id)
            sequence.append(seq)

        if not ids:
            logger.error("Could not load ids for PDB: %s, %s", pdb, chain)
            raise core.InvalidState("Could not find ids")

        if not sequence:
            logger.error("Could not load sequence for PDB: %s, %s", pdb, chain)
            raise core.InvalidState("Could not find sequence")

        sequence = ''.join(sequence)
        return {'ids': ids, 'sequence': sequence}

    def correlate(self, correlation_id, ref, target):
        if not ref or not ref.get('ids') or not ref.get('sequence'):
            logger.error("Not given complete reference: %s", ref)
            raise core.InvalidState("Incomplete reference")

        if not target or not target.get('ids') or not target.get('sequence'):
            logger.error("Not given complete target: %s", target)
            raise core.InvalidState("Incomplete target")

        data = {
            'reference': ref['sequence'],
            'reference_ids': ref['ids'],
            'target': target['sequence'],
            'target_ids': target['ids']
        }

        helper = utils.WebRequestHelper(method='post', parser=Parser())
        helper.parser.additional = {'correspondence_id': correlation_id}
        result = helper(URL, data=data, headers={'accept': 'application/json'})
        logger.info("Found %s correlations", len(result))
        return [Corr(**d) for d in result]

    def possible(self, pdb):
        return self.nr_util.members(pdb)

    def known(self, pdb):
        possible = self.possible(pdb)
        with self.session() as session:
            query = session.query(Info).\
                join(Corr, Corr.correspondence_id == Info.id).\
                filter(Info.pdb2.in_(possible))
            return set([result.pdb2 for result in query])

    def missing(self, pdb):
        return set(self.possible(pdb)) - set(self.known(pdb))

    def missing_ids(self, pdb):
        missing = self.missing(pdb)
        with self.session() as session:
            query = session.query(Info).\
                filter(Info.pdb1 == pdb).\
                filter(Info.pdb2.in_(missing))
            return [(result.id, result.pdb2) for result in query]

    def remove(self, pdb):
        possible = self.possible(pdb)
        with self.session() as session:
            query = session.query(Info.id).\
                filter(Info.pdb1 == pdb).\
                filter(Info.pdb2.in_(possible))
            ids = [id[0] for id in query]

        with self.session() as session:
            session.query(Corr).\
                filter(Corr.correspondence_id.in_(ids)).\
                delete(synchronize_session=False)

    def has_data(self, pdb):
        return not bool(self.missing(pdb))

    def data(self, pdb, **kwargs):
        data = []
        known = {}
        for corr_id, other in self.missing_ids(pdb):
            other_chain = self.st_util.longest_chain(other)
            logger.info("Using chain %s in %s", other_chain, other)

            pdb_chain = self.st_util.longest_chain(pdb)
            logger.info("Using chain %s in %s", pdb_chain, pdb)

            if pdb_chain not in known:
                known[pdb_chain] = self.structure_data(pdb_chain, pdb)
            other_data = self.structure_data(other_chain, other)
            data.extend(self.correlate(corr_id, known[pdb_chain], other_data))

        return data
