import csv
import logging
import itertools as it
import collections as coll

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import NtNtCorrespondences
from models import PdbCorrespondences
from models import PdbCoordinates
from utils import DatabaseHelper
from utils import WebRequestHelper

URL = 'http://localhost:8080/services/correlations'

CURRENT_REP_QUERY = '''
select rep_id
from nr_current_representative
where
    pdb_id = :val
'''

POLYMER_UNITS_QUERY = '''
select *
from polymer_units
where
    pdb_id = :pdb
    and model = 1
    and chain = :chain
order by polymer_id
'''


class CorrelationResponseParser(object):
    def __call__(self, text):
        reader = csv.DictReader(text)
        return [(row['reference'], row['target']) for row in reader]


class Loader(MotifAtlasBaseClass, DatabaseHelper):
    request = WebRequestHelper(method='post')
    parse = CorrelationResponseParser()

    def __init__(self, maker):
        MotifAtlasBaseClass.__init__(self)
        DatabaseHelper.__init__(self, maker)

    def reference(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_REP_QUERY, {'val': pdb})
            return set([rep[0] for rep in result.fetchall()])

    def structure_data(self, chain, pdb):
        with self.session() as session:
            result = session.execute(POLYMER_UNITS_QUERY,
                                     {'pdb': pdb, 'chain': chain})

        Record = coll.namedtuple('Record', result.keys())
        ids = []
        sequence = []
        results = it.imap(lambda r: Record(*r), result.fetchall())
        for _, nts in it.groupby(results, lambda r: r.polymer_id):
            nts = list(nts)
            ids.extend([nt.id for nt in nts])
            sequence.append(''.join([nt.unit for nt in nts]))
        return {'ids': ids, 'sequence': sequence}

    def correlate(self, correlation_id, reference, target):
        payload = {
            'reference': reference['sequence'][0],
            'reference_ids': reference['ids'],
            'target': target['sequence'],
            'target_ids': target['ids']
        }
        response = self.request(URL, payload=payload)
        parsed = self.parse(response)
        data = []
        for (ref, tar) in parsed:
            data.append(NtNtCorrespondences(unit1_id=ref,
                                            unit2_id=tar,
                                            correspondence_id=correlation_id))
        return data

    def has_correspondence(self, reference, pdb):
        with self.session() as session:
            query = session.query(PdbCorrespondences).\
                filter_by(pdb1=reference, pdb2=pdb)
            return bool(query.count())

    def correlation_id(self, reference, pdb):
        logging.debug("Adding entry for %s-%s", reference, pdb)
        data = PdbCorrespondences(pdb1=reference, pdb2=pdb)
        self.store(data)
        logging.debug("Using ID %s", data.id)
        return data.id

    def longest_chain(self, pdb):
        with self.session() as session:
            query = session.query(PdbCoordinates).\
                filter_by(pdb=pdb, model=1).\
                filter(PdbCoordinates.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(PdbCoordinates.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def data(self, pdb):
        references = self.reference(pdb)
        for reference in references:
            if not self.has_correspondence(reference, pdb):
                logging.debug("Skipping correlating with: %s", reference)
                continue

            logging.debug("Using reference: %s", reference)
            corelation_id = self.correlation_id(reference, pdb)

            ref_chain = self.longest_chain(reference)
            logging.info("Using chain %s in reference %s", ref_chain,
                         reference)

            pdb_chain = self.longest_chain(pdb)
            logging.info("Using chain %s in %s", pdb_chain, pdb)

            reference_data = self.structure_data(ref_chain, reference)
            target_data = self.structure_data(pdb_chain, pdb)

            yield self.correlate(corelation_id, reference_data,
                                 target_data)

    def __call__(self, pdbs):
        for pdb in pdbs:
            logging.info("Getting nt nt correspondence for %s", pdb)

            try:
                for data in self.data(pdb):
                    logging.info("Found %s correspondencies", len(data))
                    self.store(data)
            except:
                logging.error("Failed to store correspondencies")
