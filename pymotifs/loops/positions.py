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

class Loader(core.Loader):
    merge_data = True
    dependencies = set([UnitInfoLoader, InfoLoader])
    allow_no_data = True

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.precomputed = self.config['locations']['loops_mat_files']

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                join(mod.LoopPositions, 
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                filter(or_(mod.LoopInfo.pdb_id.in_(pdbs), mod.LoopInfo.type=='NA')).\
                distinct()
            dn_process = [r.pdb_id for r in query]
        
        to_use = sorted(set(pdbs).difference(dn_process))
            
        if not to_use
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
        """Determine the known
        """

        mapping = {}
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.PdbInfo,
                     mod.PdbInfo.pdb_id == mod.LoopInfo.pdb_id).\   
                join(mod.LoopInfo,
                     mod.LoopInfo.loop_id == mod.LoopPositions.loop_id).\
                filter(mod.PdbInfo.loops_checked==0).\  
                filter(mod.LoopInfo.pdb_id == pdb)
            
            for result in query:
                data = utils.row2dict(result)
                mapping[(result.loop_id, result.position)] = data

        return mapping

    def annotations(self, pdb, remove=True):
        """Call matlab and parse the annotations to create a list of unit id to
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

        [output_file, err_msg] = mlab.loadLoopPositions(path, nout=2)
        if err_msg != '':
            raise matlab.MatlabFailed(err_msg)

        data = self.parse(output_file)
        if remove:
            os.remove(output_file)
        return data

    def loop_units_mapping(self, pdb):
        """Load the mapping from loop id to the unit ids stored as part of the
        unit_ids field. This is used for sanity checking the given unit to loop
        positions mapping. This will load all loops mapping.

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
        """Create a callable that will normalize a comma seperated string of
        unit ids for a particular structure.

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
        """Update loop_positions table by loading data from the mat files
        stored in the PrecomputedData folder

        :param str pdb: The PDB id to get the loop positions for.
        """

        data = []
        known = self.known(pdb)

        print("loops/positions.py known")
        print(known)

        mapping = self.loop_units_mapping(pdb)

        print("loops/positions.py mapping")
        print(mapping)

        normalizer = self.normalizer(pdb)

        print("loops/positions.py mapping")
        print(normalizer)

        for row in self.annotations(pdb):
            entry = known.get((row['loop_id'], row['position']), {})
            if not entry:
                self.logger.info("New loop %s found", row['loop_id'])

            row['unit_id'] = normalizer(row['unit_id'])
            if row['loop_id'] not in mapping:
#                raise core.InvalidState("Unknown loop: %s" % row['loop_id'])
                self.logger.error("Unknown loop: %s" % row['loop_id'])
            else:
                if row['unit_id'] not in mapping[row['loop_id']]:
                    msg = "Unit %s not part of expected units %s for %s"
                    entry = (row['unit_id'], mapping[row['loop_id']], row['loop_id'])

#                    raise core.InvalidState(msg % entry)
                    self.logger.error(msg % entry)
                else:
                    entry.update(row)
                    entry['loop_id'] = row['loop_id']
                    data.append(mod.LoopPositions(**entry))

        return data
