"""Load the loop positions into the database. This will run the matlab code to
extract the loops and place them into the database.
"""

import os
import csv

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod
from pymotifs.utils import matlab
from pymotifs.utils.correct_units import Correcter
from pymotifs.units.info import Loader as UnitInfoLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from sqlalchemy import or_

class Loader(core.Loader):
    merge_data = True
    dependencies = set([UnitInfoLoader, InfoLoader])
    allow_no_data = True

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.precomputed = self.config['locations']['loops_mat_files']

    def to_process(self, pdbs, **kwargs):
        # make a list of pdb ids that don't yet have loop information in both
        # loop_info and loop_positions table
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                join(mod.LoopPositions,
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                distinct()
            # pdb ids that have entries in both loop_info and loop_positions table
            # Note: does not check that *every* loop in loop_info is also in loop_positions
            dn_process = [r.pdb_id for r in query]

        # pdb ids that don't aleady have data in both loop_info and loop_positions
        to_use = sorted(set(pdbs).difference(dn_process)) #Remove pdbs with entries in loop_positions

        # find pdb ids that have no loops at all
        # those are indicated by the presence of a 000 loop with type NA
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                filter(mod.LoopInfo.type == 'NA').\
                distinct()
            dn_process = [r.pdb_id for r in query] #list of pdbs with corresponding entries in loop_info and type='NA'

        # problemtic pdbs for now, they are wasting so much time and do nothing 2023-10-04 
        problemtic = set('1A34', '1BYX', '1ELH', '1G3A', '1H1K', '1MHK', '1N35', '1N38', '1TFW', '1TFY', '1UON', '2BGG', '2E9Z', '2F8S', '2F8T', '2G91', '2I91', '2M1O', '2M23', '2NUG', '2R92', '2RFK', '2VAL', '354D', '377D', '3AVU', '3AVV', '3AVW', '3AVX', '3AVY', '3BSN', '3BSO', '3H5X', '3H5Y', '3HTX', '3KNA', '3LRN', '3MQK', '3NCU', '3NMA', '3PLA', '3VNV', '3ZC0', '4E78', '4GV9', '4K4U', '4K4V', '4KTG', '4OQ8', '4W5O', '4W5Q', '4W5R', '4W5T')
        
        # Remove pdbs with no loops
        to_use = sorted(set(to_use).difference(dn_process))

        # Remove problemtic pbds
        to_use = sorted(set(to_use).difference(problemtic))

        if not to_use:
            raise core.Skip("Nothing to process")

        return to_use

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(mod.LoopInfo).filter_by(pdb_id=pdb)
            ids = [result.loop_id for result in query]

        if not ids:
            return True

        with self.session() as session:
            return session.query(mod.LoopPositions).\
                filter(mod.LoopPositions.loop_id.in_(ids)).\
                delete(synchronize_session=False)

        return True

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.LoopInfo,
                     mod.LoopInfo.loop_id == mod.LoopPositions.loop_id).\
                filter(mod.LoopInfo.pdb_id == pdb)
            return bool(query.count())

    def parse(self, filename):
        keys = ['loop_id', 'position', 'unit_id', 'bulge', 'flanking',
                'border']
        int_keys = ['position', 'bulge', 'flanking', 'border']
        data = []
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for row in reader:
                entry = dict(zip(keys, row))
                for key in int_keys:
                    entry[key] = int(entry[key])
                data.append(entry)
            return data

    def known(self, pdb):
        """
        Query loop_positions and return a dictionary that maps
        (loop_id,position) tuple to all the data on that row of the table
        """

        mapping = {}
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.LoopInfo,
                     mod.LoopInfo.loop_id == mod.LoopPositions.loop_id).\
                filter(mod.LoopInfo.pdb_id == pdb)

            for result in query:
                data = utils.row2dict(result)
                mapping[(result.loop_id, result.position)] = data

        return mapping

    def annotations(self, pdb, remove=True):
        """
        Call matlab and parse the annotations to create a list of unit id to
        loop mappings.

        :param str pdb: The pdb id to use.
        :param Bool remove: Flag to indicate if the produced file should be
        removed.
        :returns: The annotations produced by matlab.
        """

        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        path = str(os.path.join(self.precomputed, pdb))
        try:
            if not os.path.exists(path):
                os.mkdir(path)
        except:
            raise core.InvalidState("Could not create %s for matlab" % path)

        # run Matlab to get detailed data about each loop, save in output_file
        [output_file, err_msg] = mlab.loadLoopPositions(path, nout=2)
        if err_msg != '':
            raise matlab.MatlabFailed(err_msg)

        # read output_file and parse
        data = self.parse(output_file)

        if remove:
            # clean up the temporary output file
            os.remove(output_file)

        return data

    def loop_units_mapping(self, pdb):
        """
        Load the mapping from loop id to the unit ids stored as part of the
        unit_ids field. This is used for sanity checking the given unit to loop
        positions mapping. This will load all loops mapping.

        Queries loop_info table, unit_ids column
        Returns a dictionary mapping loop_id to the set of unit ids

        :param str pdb: The PDB id to use.
        :returns: A dictonary mapping from loop id to unit ids.
        """

        with self.session() as session:
            query = session.query(mod.LoopInfo).\
                filter_by(pdb_id=pdb)

            mapping = {}
            for result in query:
                mapping[result.loop_id] = set(result.unit_ids.split(','))
            return mapping

    def normalizer(self, pdb):
        """
        Create a callable that will normalize a comma seperated string of
        unit ids for a particular structure.

        Guess: turn modified nucleotide sequence into their parent sequence

        :param str pdb: The PDB id to use for normalization.
        :returns: A callable that will translate a comma seperated string of
        unit ids in the structure.
        """

        correcter = Correcter(self.config, self.session)
        mapping = correcter.normalized_mapping(pdb)
        def fn(unit_id):
            unit = correcter.as_unit(unit_id)
            return correcter.correct(mapping, unit)
        return fn

    def data(self, pdb, **kwargs):
        """
        Add to loop_positions table by loading data from the mat files
        stored in the PrecomputedData folder

        Note:  Matlab is not able to read modified nucleotides

        :param str pdb: The PDB id to get the loop positions for.
        """

        data = []

        # Query loop_positions and return a dictionary that maps
        # (loop_id,position) tuple to all the data on that row of the table
        known = self.known(pdb)

        # mapping is a dictionary mapping loop_id to the set of unit ids
        mapping = self.loop_units_mapping(pdb)

        # Guess:  turn modified nucleotide sequence into parent sequence
        normalizer = self.normalizer(pdb)

        # loop over
        for row in self.annotations(pdb):
            entry = known.get((row['loop_id'], row['position']), {})
            if not entry:
                self.logger.info("New loop %s found", row['loop_id'])

            row['unit_id'] = normalizer(row['unit_id'])
            if row['loop_id'] not in mapping:
                self.logger.error("Unknown loop: %s" % row['loop_id'])
            else:
                if row['unit_id'] not in mapping[row['loop_id']]:
                    msg = "Unit %s not part of expected units %s for %s"
                    entry = (row['unit_id'], mapping[row['loop_id']], row['loop_id'])

                    self.logger.error(msg % entry)
                else:
                    entry.update(row)
                    entry['loop_id'] = row['loop_id']
                    data.append(mod.LoopPositions(**entry))

        return data
