"""
This is a module to extract loops from structures. It uses matlab to find all
loops and then will save them in the correct location as specificed by
'locations'. It also stores the the loop information into the database.
"""

from pymotifs.models import LoopsAll
from pymotifs import core
from pymotifs import utils


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

    def _extract_loops(self, pdb_id, loop_type):
        """Loops - array of FR3D File structures. l - its length"""

        [Loops, l, err_msg] = self.matlab.extractLoops(pdb_id, loop_type,
                                                       nout=3)

        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        if Loops == 0:
            raise core.SkipValue('No %s in %s' % loop_type, pdb_id)

        self.logger.info('Found %i loops', l)
        return (Loops, l)

    def _get_loop_mapping(self, pdb_id, loop_type):
        """Compute a mapping from the nts to the loop id.  This is used
        for setting ids by either looking up the old known id or creating a new
        one if no one is found.

        :pdb_id: The pdb id to search.
        :loop_type: The loop type.
        :returns: A dictonary mapping from the nts to the loop id. If there are
        no loops it returns an empty dictonary.
        """

        mapping = {}
        with self.session() as session:
            query = self.query(session, (pdb_id, loop_type))
            for result in query:
                mapping[result.nt_ids] = utils.row2dict(result)
        return mapping

    def _next_loop_number_string(self, current):
        """Compute the next loop number string. This will pad to either 3 or 6
        characters with zeros. If the next number is over 999 we use 6,
        otherwise 3 as the length to pad to.

        :current: The current loop count.
        :returns: A string of the next loop id.
        """

        next_number = current + 1
        if next_number > 999:
            return str(next_number).rjust(6, '0')
        return str(next_number).rjust(3, '0')

    def _get_loop_id(self, nts, pdb_id, loop_type, mapping, count):
        """Compute the loop id to use for the given unit string. This will
        build a string like IL_1S72_001 or IL_4V4Q_001000. In structures with
        over 999 loops, we will pad with zeros to 6 characters, but keep the
        stanadrd padding to 3 characters otherwise.

        :nts: The concanated unit string.
        :pdb_id: The pdb id to use.
        :loop_type: The type of loop.
        :mapping: A mapping from unit string to known loop_id.
        :count: The current number of known loops.
        :returns: A string of the new loop id.
        """

        if nts in mapping:
            self.logger.debug('Nucleotides %s matched %s', nts, mapping[nts])
            return mapping[nts]

        # format examples: IL_1S72_001, IL_4V4Q_001000
        str_number = self._next_loop_number_string(count)
        loop_id = '_'.join([loop_type, pdb_id, str_number])
        self.logger.info('Created new loop id %s, for nucleotides %s',
                         loop_id, nts)
        return loop_id

    def _loop_objects(self, loops, l, pdb_id, loop_type, mapping):
        data = []
        for i in xrange(l):
            loop = loops[i].LoopsAll_table
            loop_id = self._get_loop_id(loop.full_id, pdb_id, mapping,
                                        len(data))

            loops[i].Filename = loop_id
            data.append(LoopsAll(
                id=loop_id,
                type=loop_type,
                pdb=pdb_id,
                sequential_id=loop_id.split("_")[-1],
                length=int(loops[i].NumNT[0][0]),
                seq=loop.seq,
                r_seq=loop.r_seq,
                nwc_seq=loop.nwc,
                r_nwc_seq=loop.r_nwc,
                nt_ids=loop.full_id,
                loop_name=loop.loop_name))

        return data

    def _save_mat_files(self, loops):
        """Pass the Loops structure array back to matlab so that it can
        save the .mat files in the specified location.
        """

        location = self.config['locations']['loops_mat_files']
        [status, err_msg] = self.matlab.aSaveLoops(loops, location, nout=2)

        if status != 0:
            raise core.StageFailed("Could not save all loop mat files")

        self.logger.info('mat files saved')

    def data(self, entry, **kwargs):
        pdb_id, loop_type = entry
        (loops, l) = self._extract_loops(pdb_id, loop_type)
        self._save_mat_files(loops)
        mapping = self._get_loop_mapping(self, pdb_id, loop_type)
        return self._loop_objects(loops, l, pdb_id, loop_type, mapping)
