"""
This is a module to extract loops from structures. It uses matlab to find all
loops and then will save them in the correct location as specificed by
'locations'. It also stores the the loop information into the database.
"""
import os

from pymotifs.models import LoopInfo
from pymotifs import core
from pymotifs.utils.units import Translator
from pymotifs.loops.release import Loader as ReleaseLoader


class Loader(core.SimpleLoader):
    loop_types = ['IL', 'HL', 'J3']
    merge_data = True
    allow_no_data = True
    dependencies = set([ReleaseLoader])

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)

    def must_recompute(self, entry, recalculate=False, **kwargs):
        return bool(recalculate) or \
            self.config[self.name].get('recompute') or \
            self.config['recalculate'].get(entry[1])

    def query(self, session, pdb):
        return session.query(LoopInfo).filter_by(pdb_id=pdb)

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

    def _get_loop_id(self, units, pdb_id, loop_type, mapping):
        """Compute the loop id to use for the given unit string. This will
        build a string like IL_1S72_001 or IL_4V4Q_001000. In structures with
        over 999 loops, we will pad with zeros to 6 characters, but keep the
        stanadrd padding to 3 characters otherwise.

        :units: The concanated unit or nt id string.
        :pdb_id: The pdb id to use.
        :loop_type: The type of loop.
        :mapping: A mapping from unit string to known loop_id.
        :count: The current number of known loops.
        :returns: A string of the new loop id.
        """

        if units in mapping:
            self.logger.debug('Nucleotides %s matched %s', units,
                              mapping[units])
            return mapping[units]

        # format examples: IL_1S72_001, IL_4V4Q_001000
        str_number = self._next_loop_number_string(len(mapping))
        loop_id = '_'.join([loop_type, pdb_id, str_number])
        self.logger.info('Created new loop id %s, for nucleotides %s',
                         loop_id, units)
        return loop_id

    def _extract_loops(self, pdb_id, loop_type, mapping):
        """
        Uses matlab to extract the loops for a given structure of a specific
        type. This will also save the loop files into the correct place.
        """

        location = os.path.join(self.config['locations']['loops_mat_files'])
        try:
            mlab = core.Matlab(self.config['locations']['fr3d_root'])
            [loops, count, err_msg] = \
                mlab.extractLoops(pdb_id, loop_type, nout=3)
        except Exception as err:
            self.logger.exception(err)
            raise err

        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        if loops == 0:
            self.logger.warning('No %s in %s', loop_type, pdb_id)
            return []

        self.logger.info('Found %i loops', count)

        data = []
        for index in xrange(count):
            loop = loops[index].AllLoops_table
            loop_id = self._get_loop_id(loop.loop_name, pdb_id, loop_type,
                                        mapping)
            mapping[loop.full_id] = str(loop_id)
            loops[index].Filename = str(loop_id)

            data.append(LoopInfo(
                loop_id=loop_id,
                type=loop_type,
                pdb=str(pdb_id),
                sequential_id=loop_id.split("_")[-1],
                length=int(loops[index].NumNT[0][0]),
                seq=loop.seq,
                r_seq=loop.r_seq,
                nwc_seq=loop.nwc,
                r_nwc_seq=loop.r_nwc,
                unit_ids=loop.full_id,
                loop_name=loop.loop_name))
        try:
            [status, err_msg] = mlab.aSaveLoops(loops, str(location), nout=2)
        except Exception as err:
            self.logger.exception(err)
            raise err

        if status != 0:
            raise core.MatlabFailed("Could not save all loop mat files")

        self.logger.info("Saved %s_%s loop mat files", loop_type, pdb_id)

        return data

    def _mapping(self, pdb_id, loop_type):
        """Compute a mapping from the nts to the loop id.  This is used
        for setting ids by either looking up the old known id or creating a new
        one if no one is found.

        :pdb_id: The pdb id to search.
        :loop_type: The loop type.
        :returns: A dictonary mapping from the nts to the loop id. If there are
        no loops it returns an empty dictonary.
        """

        mapping = {}
        translator = Translator(self.session.maker)

        with self.session() as session:
            query = self.query(session, pdb_id)
            for result in query:
                unit_list = result.nt_ids.split(',')
                seperator = unit_list[0][4]

                units = None
                if seperator == '_':
                    self.logger.debug("Translating %s to unit ids",
                                      result.nt_ids)
                    unit_list = translator.translate(unit_list)
                    units = ','.join(unit_list)
                    self.logger.debug("Translated to: %s", units)
                elif seperator == '|':
                    self.logger.debug("No need to translate unit ids")
                    units = result.nt_ids
                else:
                    raise ValueError("Unknown seperator type: %s", seperator)

                mapping[units] = result.id
        return mapping

    def data(self, pdb, **kwargs):
        data = []
        for loop_type in self.loop_types:
            mapping = self._mapping(pdb, loop_type)
            data.extend(self._extract_loops(pdb, loop_type, mapping))

        if not data:
            raise core.Skip("No loops found in %s", pdb)

        return data
