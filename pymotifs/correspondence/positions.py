"""
This is a stage to align two experimental sequences and store the alignment
between each nucleotide.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.correspondence.info import Loader as CorrLoader
from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader

from pymotifs.utils.alignment import align
from pymotifs.utils.alignment import one_to_one_alignment


class Loader(core.Loader):
    """
    A loader for computing the position to position alignment and storing
    it.
    """

    mark = False
    dependencies = set([CorrLoader, InfoLoader, PositionLoader])

    def to_process(self, pdbs, **kwargs):
        """
        Ignore the list of pdbs passed in.
        Query the correspondence_info table.
        Get all correspondence_id values without a length,
        indicating that no correspondence is calculated yet.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo.correspondence_id).\
                filter(mod.CorrespondenceInfo.length == None)
            if not query.count():
                raise core.Skip("Skipping positions, no new correspondences")
            return [result.correspondence_id for result in query]


    def has_data(self, corr_id, **kwargs):
        """Check if we have data for the given correspondence id. This will
        check if the stored length is not null, in which case it has been
        computed and summarized before. This can cause problems if this stage
        is rerun before doing a summary and summary.

        :param int corr_id: The correspondence id to use.
        :returns: A boolean.
        """

        with self.session() as session:
            corr = session.query(mod.CorrespondenceInfo).\
                filter_by(correspondence_id=corr_id)

            if not corr.count():
                return False

            corr = corr.one()
            return (corr.length != None) and (corr.length != 0)

    def remove(self, corr_id, **kwargs):
        """Remove all the data for the given correspondence id.

        :param int corr_id: The correspondence id to cleanup.
        """

        if kwargs.get('dry_run'):
            return True

        with self.session() as session:
            session.query(mod.CorrespondencePositions).\
                filter_by(correspondence_id=corr_id).\
                delete(synchronize_session=False)

    def sequence(self, exp_id):
        """Load all information about the experimental sequence with the given
        id. This will load both the ids and the sequence. The ids will be a
        list of numbers, while the sequence is a string.

        :param int exp_id: The experimental sequence id.
        :returns: A dictionary of the ids and sequence for the given id.
        """

        ids = []
        sequence = []
        with self.session() as session:
            query = session.query(mod.ExpSeqPosition).\
                filter(mod.ExpSeqPosition.exp_seq_id == exp_id).\
                order_by(mod.ExpSeqPosition.index)

            if not query.count():
                raise core.InvalidState("Could not get sequence for %s" % exp_id)

            for index, result in enumerate(query):
                seq_id = result.exp_seq_position_id
                seq = result.normalized_unit or 'N'
                ids.append(seq_id)
                sequence.append(seq)

        return {'ids': ids, 'sequence': ''.join(sequence)}

    def info(self, corr_id):
        """Look up the sequences used in some correspondence.

        :param int corr_id: The id of the correspondence to lookup.
        :returns: A tuple of experimental sequences used.
        """

        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo).\
                filter_by(correspondence_id=corr_id)
            if not query.count():
                raise core.InvalidState("Unknown correspondence %s" % corr_id)

            result = query.one()
            return (result.exp_seq_id_1, result.exp_seq_id_2)

    def align_sequences(self, corr_id, ref, target): ## rename !! align_sequences
        """Run the alignment on two sequences. This will do an alignment are
        return lists that contain the ids for aligned positions only.

        :param int corr_id: The correspondence id to use.
        :param dict ref: The reference sequence.
        :param dict target: The target sequence.
        :returns: A list of dictionaries for each position in the alignment. It
        lists which positions are aligned.
        """

        seq1, seq2 = self.info(corr_id)
        self.logger.info('show seq1: %d and seq2: %d'%(seq1,seq2))
        self.logger.info('show the corr_id: %d'%corr_id)
        with self.session() as session:
            query = session.query(mod.ExpSeqInfo.entity_type).\
                filter(mod.ExpSeqInfo.exp_seq_id.in_([seq1,seq2]))
        # for result in query:
        #     self.logger.info("show the query result:%s %s %s" % (result.exp_seq_id, result.entity_type, corr_id))
        #     print(result.entity_type)
        #     print(result.exp_seq_id)


        entity_type_check = set()
        # this is a double check for entity types because we have checked sequence pairs when we are making sequence pairs.
        for result in query:
            
            if result.entity_type == 'rna':
                entity_type_check.add('rna')
            elif result.entity_type == 'dna':
                entity_type_check.add('dna')
            elif result.entity_type == 'hybrid':
                entity_type_check.add('hybrid')

        
        if len(entity_type_check) > 1:
            raise core.InvalidState('The entity types of the sequence pair are not identical')
        elif list(entity_type_check) == ['rna']:
            results = align([ref, target])

            data = []
            for index, result in enumerate(results):

                data.append({
                    'exp_seq_position_id_1': result[0],
                    'exp_seq_position_id_2': result[1],
                    'correspondence_id': corr_id,
                    'index': index
                })
            return data
            
        elif list(entity_type_check) == ['dna']:
            results = one_to_one_alignment([ref, target])
            data = []
            for index, result in enumerate(results):
                data.append({
                    'exp_seq_position_id_1': result[0],
                    'exp_seq_position_id_2': result[1],
                    'correspondence_id': corr_id,
                    'index': index
                })
            return data
           
        elif list(entity_type_check) == ['hybrid']:
            raise core.InvalidState('The hybrid pair alignments have not defined.')
        
        # return data


    # def info(self, corr_id):
    #     """Look up the sequences used in some correspondence.

    #     :param int corr_id: The id of the correspondence to lookup.
    #     :returns: A tuple of experimental sequences used.
    #     """

    #     with self.session() as session:
    #         query = session.query(mod.CorrespondenceInfo).\
    #             filter_by(correspondence_id=corr_id)
    #         if not query.count():
    #             raise core.InvalidState("Unknown correspondence %s" % corr_id)

    #         result = query.one()
    #         return (result.exp_seq_id_1, result.exp_seq_id_2)

    def align_sequences_old(self, corr_id, ref, target): ## rename !! align_sequences
        """Run the alignment on two sequences. This will do an alignment are
        return lists that contain the ids for aligned positions only.

        :param int corr_id: The correspondence id to use.
        :param dict ref: The reference sequence.
        :param dict target: The target sequence.
        :returns: A list of dictionaries for each position in the alignment. It
        lists which positions are aligned.
        """

        results = align([ref, target])

        data = []
        for index, result in enumerate(results):

            data.append({
                'exp_seq_position_id_1': result[0],
                'exp_seq_position_id_2': result[1],
                'correspondence_id': corr_id,
                'index': index
            })
        return data



    def data(self, corr_id, **kwargs):
        """Compute the positions for a correspondence. This will do a position
        to position correspondence for all

        :param int corr_id: The correspondence id.
        :yields: The correspondences by positions in both directions.
        """

        exp_id1, exp_id2 = self.info(corr_id)
        sequence1 = self.sequence(exp_id1)
        sequence2 = self.sequence(exp_id2)
        ## self.logger.info('corr_id: %s  !!!!!!!!!!!!!%s!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'%(corr_id,sequence1))
        alignment = self.align_sequences(corr_id, sequence1, sequence2)      ###
        ## self.logger.info('alignment %s'%alignment)
        for position in alignment:
            yield mod.CorrespondencePositions(**position)

            pos1 = position['exp_seq_position_id_1']
            pos2 = position['exp_seq_position_id_2']
            if pos1 != pos2:
                rev = dict(position)
                rev['exp_seq_position_id_1'] = pos2
                rev['exp_seq_position_id_2'] = pos1

                yield mod.CorrespondencePositions(**rev)
