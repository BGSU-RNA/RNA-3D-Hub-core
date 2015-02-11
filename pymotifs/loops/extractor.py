"""
This is a module to extract loops from structures. It uses matlab to find all
loops and then will save them in the correct location as specificed by
'locations'

"""

from pymotifs.models import LoopsAll
from pymotifs import core


class Loader(core.SimpleLoader):
    loop_types = ['IL', 'HL', 'J3']
    merge_data = True
    allow_no_data = True

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.matlab = core.Matlab(self.config['locations']['fr3d_root'])

    def transform(self, pdb, **kwargs):
        return [(pdb, type) for type in self.loop_types]

    def must_recompute(self, entry, recalculate=False, **kwargs):
        return bool(recalculate) or \
            self.config[self.name].get('recompute') or \
            self.config['recalculate'].get(entry[1])

    def query(self, session, entry):
        return session.query(LoopsAll).\
            filter_by(pdb=entry[0], type=entry[1])

    def data(self, entry, **kwargs):
        pdb_id, loop_type = entry
        (loops, l) = self.extract_loops(pdb_id, loop_type)
        # self.save_mat_files(loops)
        return self.loop_objects(loops, l, pdb_id, loop_type)

    def extract_loops(self, pdb_id, loop_type):
        """Loops - array of FR3D File structures. l - its length"""

        [Loops, l, err_msg] = self.matlab.extractLoops(pdb_id, loop_type,
                                                       nout=3)

        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        if Loops == 0:
            self.logger.info('No %s in %s', loop_type, pdb_id)
            return (0, 0)
        else:
            self.logger.info('Found %i loops', l)
            return (Loops, l)

    def loop_objects(self, loops, l, pdb_id, loop_type):
        data = []
        for i in xrange(l):
            loop_id = self._get_loop_id(loops[i].AllLoops_table.full_id,
                                        pdb_id, loop_type)
            loops[i].Filename = loop_id
            data.append(LoopsAll(
                id=loop_id,
                type=loop_type,
                pdb=pdb_id,
                sequential_id=loop_id[-3:],
                length=int(loops[i].NumNT[0][0]),
                seq=loops[i].AllLoops_table.seq,
                r_seq=loops[i].AllLoops_table.r_seq,
                nwc_seq=loops[i].AllLoops_table.nwc,
                r_nwc_seq=loops[i].AllLoops_table.r_nwc,
                nt_ids=loops[i].AllLoops_table.full_id,
                loop_name=loops[i].AllLoops_table.loop_name))
        return data

    def save_mat_files(self, loops):
        """Pass the Loops structure array back to matlab so that it can
        save the .mat files in the specified location.
        """

        location = self.config['locations']['loops_mat_files']
        [status, err_msg] = self.matlab.aSaveLoops(loops, location, nout=2)

        if status == 0:
            self.logger.info('mat files saved')
        else:
            raise core.StageFailed("Could not save all loop mat files")

    def next_loop_number_string(self, current):
        next_number = current + 1
        if next_number > 999:
            return str(next_number).rjust(6, '0')
        return str(next_number).rjust(3, '0')

    def _get_loop_id(self, full_id, pdb_id, loop_type):
        """Get a loop id for the loop. If its already known we us the current
        ID, otherwise we create a new one.

        :full_id: A string of the list of nucleotides in 5'->3' order joined by
        commas.
        :pdb_id: The pdb id we are getting a loop id for.
        :loop_type: The loop type, such as IL, HL, or J3.
        :returns: A string representing the loop id, may be new.
        """

        with self.session() as session:
            loop = session.query(LoopsAll).\
                filter(LoopsAll.nt_ids == full_id).\
                first()

            if loop:
                self.logger.debug('Full_id %s matched %s', full_id, loop.id)
                return loop.id

        # count the loops already in the db
        with self.session() as session:
            seq_id = session.query(LoopsAll).\
                filter(LoopsAll.pdb == pdb_id).\
                filter(LoopsAll.type == loop_type).\
                count()

            # format example: IL_1S72_001
            str_number = self.next_loop_number_string(seq_id)
            id = '_'.join([loop_type, pdb_id, str_number])
            self.logger.info('Created new loop id %s', id)
            return id


# SELECT * FROM loops_all AS t1
# LEFT JOIN motifversions.`loops_all` AS t2
# ON t1.id=t2.id
# WHERE t1.type!=t2.type
# OR t2.pdb!=t1.pdb
# OR t1.`sequential_id`!=t2.`sequential_id`
# OR t1.length!=t2.length
# OR t1.seq!=t2.seq
# OR t1.`r_nwc_seq`!=t2.`r_nwc_seq`
# OR t1.`r_seq`!=t2.`r_seq`
# OR t1.`nt_ids`!=t2.`nt_ids`
# OR t1.`loop_name` !=t2.`loop_name`
# OR t1.`pdb_file`!=t2.`pdb_file`
# OR t1.`nwc_seq`!=t2.`nwc_seq` #junctions are different;
