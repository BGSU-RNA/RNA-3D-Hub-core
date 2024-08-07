"""This is a module to extract loops from structures. It uses matlab to find
all loops and then will save them in the correct location as specificed by
'locations'. It also stores the the loop information into the database.
"""

import os
import functools as ft

from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod
from pymotifs.utils.correct_units import Correcter

from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.mat_files import Loader as MatLoader
from pymotifs.interactions.loader import Loader as InteractionLoader

class Loader(core.SimpleLoader):
    loop_types = ['IL', 'HL', 'J3']
    merge_data = True
    allow_no_data = True
    dependencies = set([PdbLoader, MatLoader, InteractionLoader, UnitLoader])
    save_loops = True

    def to_process(self, pdbs, **kwargs):

        # prevent trying to extract loops with Matlab when filling in DNA releases
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        self.logger.info("nr_molecule_parent_current: %s" % nr_molecule_parent_current)
        if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
            raise core.Skip("loops.extractor does not annotate DNA structures")

        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                distinct()
            known = [r.pdb_id for r in query]

        to_use = sorted(set(pdbs).difference(known))  #We want to process ONLY the pdbs that are NOT in loop_info

        self.logger.info("Extracting loops from %s" % to_use)

        if not to_use:
            raise core.Skip("no new PDB ids that need loops extracted")
        return to_use

    def query(self, session, pdb):
        return session.query(mod.LoopInfo).filter_by(pdb_id=pdb)

    def remove(self, *args, **kwargs):
        """Does not actually remove from the DB. We always want the loop ids to
        be consistent so we do not automatically remove loops.
        """

        self.logger.info("We don't actually remove data for loop extractor")

    def normalizer(self, pdb):
        """Create a callable that will normalize a comma seperated string of
        unit ids for a particular structure.

        :param str pdb: The PDB id to use for normalization.
        :returns: A callable that will translate a comma seperated string of
        unit ids in the structure.
        """

        corrector = Correcter(self.config, self.session)
        mapping = corrector.normalized_mapping(pdb)

        def fn(unit_string):
            unit_ids = unit_string.split(',')
            units = [corrector.as_unit(uid) for uid in unit_ids]
            corrected = corrector.correct_structure(pdb, mapping, units,
                                                    require_all=True)
            return tuple(sorted(unit for unit in corrected))

        return fn

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
        build a string like IL_1S72_001 or IL_4V4Q_001. In structures with
        over 999 loops, we will pad with zeros to 6 characters, but keep the
        stanadrd padding to 3 characters otherwise.

        :units: The concanated unit or nt id string.
        :pdb_id: The pdb id to use.
        :loop_type: The type of loop.
        :mapping: A mapping from unit string to known loop_id.
        :count: The current number of known loops.
        :returns: A string of the new loop id.
        """

        # format examples: IL_1S72_001, IL_4V4Q_001
        if units not in mapping:
            str_number = self._next_loop_number_string(len(mapping))
            loop_id = '%s_%s_%s' % (loop_type, pdb_id, str_number)
            self.logger.info('Created new loop id %s, for nucleotides %s',
                             loop_id, units)
            mapping[units] = str(loop_id)

        return mapping[units]

    def _get_fake_loop_id(self, pdb_id, loop_type):

        loop_id = '%s_%s_%s' % (loop_type, pdb_id, '000')
        self.logger.info('Created new fake loop id %s for pdb_id %s', loop_id, pdb_id)

        return loop_id

    def _extract_loops(self, pdb, loop_type, mapping, normalize):
        """Uses matlab to extract the loops for a given structure of a specific
        type. This will also save the loop files into the correct place.

        :param str pdb: PDB file to process
        :param str loop_type: The type of loop (IL, HL, J3, ...) to extract
        loops for.
        :param dict mapping: A mapping of unit ids to known loop names.
        :returns: The extracted loops.
        """
        try:
            mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
            [loops, count, err_msg] = mlab.extractLoops(pdb, loop_type, nout=3)
        except Exception as err:
            self.logger.exception(err)
            raise err

        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        if loops == 0:
            self.logger.warning('No %s in %s', loop_type, pdb)
            loop_id = self._get_fake_loop_id(pdb, loop_type)
            return [mod.LoopInfo(loop_id=loop_id,
                type = 'NA',
                pdb_id=pdb,
                sequential_id='000',
                length=0,
                seq='',
                r_seq='',
                nwc_seq='',
                r_nwc_seq='',
                unit_ids='',
                loop_name='')]

        self.logger.info('Found %i %s loops', count, loop_type)

        data = []
        for index in xrange(count):
            loop = loops[index].AllLoops_table
            full_id = normalize(loop.full_id)
            loop_id = self._get_loop_id(full_id, pdb, loop_type, mapping)
            loops[index].Filename = loop_id

            data.append(mod.LoopInfo(
                loop_id=loop_id,
                type=loop_type,
                pdb_id=pdb,
                sequential_id=loop_id.split("_")[-1],
                length=int(loops[index].NumNT[0][0]),
                seq=loop.seq,
                r_seq=loop.r_seq,
                nwc_seq=loop.nwc,
                r_nwc_seq=loop.r_nwc,
                unit_ids=','.join(full_id),
                loop_name=loop.loop_name))

        if self.save_loops:
            self.__save__(loops, self.config['locations']['loops_mat_files'])

        return data

    def __save__(self, loops, location):
        """Save the loops to a file.

        :loops: The loops matlab proxy object to save.
        :param str location: The location to saave to.
        """

        if not os.path.isdir(location):
            os.makedirs(location)

        try:
            mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
            [status, err_msg] = mlab.aSaveLoops(loops, location, nout=2)
        except Exception as err:
            self.logger.exception(err)
            raise err

        if status != 0:
            self.logger.error(mlab.last_stdout)
            raise matlab.MatlabFailed("Could not save all loop mat files")

        self.logger.debug("Saved loop mat files")

    def _mapping(self, pdb_id, loop_type, normalizer):
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
            query = self.query(session, pdb_id).filter_by(type=loop_type)
            for result in query:
                unit_ids = normalizer(result.unit_ids)
                if unit_ids in mapping:
                    self.logger.error("Loop %s duplicates %s",
                                      result.loop_id, mapping[unit_ids])
                    continue
                mapping[unit_ids] = result.loop_id
        return mapping

    def data(self, pdb, **kwargs):
        """Compute the loops in the given structure.

        :param str pdb: The structure to get loops for.
        :returns: A list of all the loops.
        """

        data = []
        normalizer = self.normalizer(pdb)
        for loop_type in self.loop_types:
            mapping = self._mapping(pdb, loop_type, normalizer)
            data.extend(self._extract_loops(pdb, loop_type, mapping, normalizer))
        return data
