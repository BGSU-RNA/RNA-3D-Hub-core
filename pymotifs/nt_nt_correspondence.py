import sys
import json
import logging
import traceback
import itertools as it
import collections as coll

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import NtNtCorrespondences as Corr
from models import PdbCorrespondences
from models import PdbCoordinates
from utils import DatabaseHelper
from utils import WebRequestHelper
from utils import EmptyResponse

logger = logging.getLogger(__name__)

URL = 'http://localhost:8080/api/services/correlations'

CURRENT_REP_QUERY = '''
select rep_id
from nr_current_representative
where
    pdb_id = :val
'''

POLYMER_UNITS_QUERY = '''
select *
from polymer_units as P
join pdb_unit_ordering as O
on
    O.nt_id = P.id
where
    P.pdb_id = :pdb
    and P.model = 1
    and P.chain = :chain
    and P.model = :model
order by P.polymer_id, O.index
'''


class CorrelationResponseParser(object):
    def __init__(self, additional=None):
        self.additional = additional or {}

    def __call__(self, text):
        if not text:
            raise EmptyResponse()

        data = []
        for row in json.loads(text):
            entry = {'unit1_id': row['unit1'], 'unit2_id': row['unit2']}
            entry.update(self.additional)
            data.append(entry)
        return data


class StructureUtil(DatabaseHelper):
    def longest_chain(self, pdb, model=1):
        with self.session() as session:
            query = session.query(PdbCoordinates).\
                filter_by(pdb=pdb, model=model).\
                filter(PdbCoordinates.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(PdbCoordinates.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def representative(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_REP_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep

    def polymer_sequences(self, pdb, chain, model=1):
        results = self.__polymer_units__(pdb=pdb, chain=chain, model=model)
        sequence = []
        for _, nts in it.groupby(results, lambda r: r.polymer_id):
            sequence.append(''.join([nt.unit for nt in nts]))
        return sequence

    def unit_ids(self, pdb, chain, model=1):
        results = self.__polymer_units__(pdb=pdb, chain=chain, model=model)
        return [result.id for result in results]

    def __polymer_units__(self, **data):
        with self.session() as session:
            result = session.execute(POLYMER_UNITS_QUERY, data)

        Record = coll.namedtuple('Record', result.keys())
        return it.imap(lambda r: Record(*r), result.fetchall())


class Loader(MotifAtlasBaseClass, DatabaseHelper):
    request = WebRequestHelper(method='post',
                               parser=CorrelationResponseParser())

    def __init__(self, maker):
        MotifAtlasBaseClass.__init__(self)
        DatabaseHelper.__init__(self, maker)
        self.util = StructureUtil(maker)

    def structure_data(self, chain, pdb):
        return {
            'ids': self.util.unit_ids(pdb, chain),
            'sequence': self.util.polymer_sequences(pdb, chain),
            'pdb': pdb
        }

    def correlation_id(self, reference, pdb):
        data = PdbCorrespondences(pdb1=reference, pdb2=pdb)
        self.store([data])
        with self.session() as session:
            data = session.query(PdbCorrespondences).\
                filter_by(pdb1=reference, pdb2=pdb).\
                one().id
        return data

    def correlate(self, reference, target):
        payload = {
            'reference': reference['sequence'][0],
            'reference_ids': reference['ids'],
            'target': target['sequence'],
            'target_ids': target['ids']
        }
        headers = {'accept': 'application/json'}
        correlation_id = self.correlation_id(reference['pdb'], target['pdb'])
        self.request.parser.additional = {'correspondence_id': correlation_id}
        response = self.request(URL, data=payload, headers=headers)
        return [Corr(**d) for d in response]

    def __base_info_query__(self, session, reference, pdb):
        # TODO: Do something to store correlation method data
        return session.query(PdbCorrespondences).\
            filter_by(pdb1=reference, pdb2=pdb)

    def has_correspondence(self, reference, pdb):
        with self.session() as session:
            query = self.__base_info_query__(session, reference, pdb)
            return bool(query.count())

    def remove_old(self, reference, pdb):
        with self.session() as session:
            query = self.__base_info_query(session, reference, pdb)
            corr_id = query.one().id
            session.query(Corr).filter_by(correspondence_id=corr_id).delete()
            query.delete()

    def data(self, pdb, recalculate=False, **kwargs):
        for reference in self.util.representative(pdb):
            if self.has_correspondence(reference, pdb):
                if not recalculate:
                    logger.debug("Skipping correlating with: %s", reference)
                    continue
                else:
                    logger.debug("Deleting old correspondences for %s, %s",
                                 reference, pdb)
                    self.remove_old(reference, pdb)

            logger.debug("Using reference: %s", reference)

            ref_chain = self.util.longest_chain(reference)
            logger.info("Using chain %s in reference %s", ref_chain, reference)

            pdb_chain = self.util.longest_chain(pdb)
            logger.info("Using chain %s in %s", pdb_chain, pdb)

            yield self.correlate(self.structure_data(ref_chain, reference),
                                 self.structure_data(pdb_chain, pdb))

    def __call__(self, pdbs, **kwargs):
        if not pdbs:
            raise Exception("No pdbs given")

        for pdb in pdbs:
            logger.info("Getting nt nt correspondence for %s", pdb)

            try:
                for data in self.data(pdb, **kwargs):
                    logger.info("Found %s correspondencies", len(data))
                    self.store(data)
            except:
                logger.error("Failed to store correspondencies for %s", pdb)
                logger.error(traceback.format_exc(sys.exc_info()))


if __name__ == '__main__':
    from utils import main
    main(Loader)
