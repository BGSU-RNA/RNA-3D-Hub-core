"""
Import all pairwise interactions.

Runs Matlab fr3d on the given files to determine all interactions and then imports
them into the database.
"""

import re
import os
import csv
import operator as op
import collections as coll

from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod

from pymotifs.mat_files import Loader as MatLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.export.cifatom import Exporter as CifAtom

IGNORE = set(['perp', 'nbif', 'bif', 'rib', 'nRib', 'rIB'])


class Loader(core.SimpleLoader):
    """A loader to generate and import the interaction annotations for
    structures.
    """

    allow_no_data = True

    dependencies = set([MatLoader, UnitLoader, PdbLoader, CifAtom])
    @property
    def table(self):
        return mod.UnitPairsInteractions

    def query(self, session, pdb):
        """Create a query to access interaction data for the given pdb.

        :session: The database session to use.
        :pdb: The pdb id
        :returns: A query to get interaction data.
        """
        ## stopped using the query function to check if we have reords for this table.
        ## we added to_process function to do the same job.
        # return 0
        # filter by program='matlab' because that is what this program is about
        return session.query(mod.UnitPairsInteractions).filter_by(pdb_id=pdb,program='matlab')

    def to_process(self, pdbs, **kwargs):
        """
        Query to find pdb ids with 'placeholder' unit ids, which means that the
        structure was processed to find pairwise interactions, none were found,
        and now we can avoid checking again.
        """

        # prevent trying to create Matlab interactions when filling in DNA releases
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        self.logger.info("nr_molecule_parent_current: %s" % nr_molecule_parent_current)

        if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
            raise core.Skip("interactions.pairwise does not annotate DNA structures")

        # find all unique pdb ids that have been processed and have a pairwise interaction or have a placeholder
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions.pdb_id).\
                filter(mod.UnitPairsInteractions.program == 'matlab').\
                distinct()
            done = set()
            for result in query:
                done.add(result.pdb_id)

        self.logger.info('Found %d pdbs with matlab pairwise interactions' % len(done))

        needed = sorted(list(set(pdbs)-done))

        self.logger.info('Found %d pdbs that need to be processed for pairwise interactions' % len(needed))

        if len(needed) == 0:
            raise core.Skip("All pdbs have been processed for pairwise interactions")

        # remove the pdb ids just found from the list of those that need to be processed
        return needed

    def interaction_type(self, family):
        """Determine the interaction type of the given interaction. This will
        return the column name in the table this should be added to. If it
        matches nothing a warning is logged and None is returned.

        :family: The interaction annotation to get the family for.
        :returns: The type of the interaction.
        """

        if re.match(r'n?s[53]{2}$', family):
            return 'f_stacks'
        elif re.match(r'n?\dBR$', family):
            return 'f_brbs'
        elif re.match(r'^n?\dBPh$', family):
            return 'f_bphs'
        elif re.match(r'^n?[ct][WHS]{2}$', family) or family == 'wat':
            return 'f_lwbp'
        elif family in IGNORE:
            return None

        self.logger.warning("Unknown interaction: %s", family)
        return None

    def parse(self, filename, pdb):
        """Reads the csv file, imports all interactions, deletes the file when
        done to avoid stale data and free up disk space

        :filename: The input filename.
        :pdb: The pdb id.
        :returns: A list of Interaction objects.
        """

        data = coll.defaultdict(dict)
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for index, row in enumerate(reader):
                if not row[0] or not row[1]:
                    msg = "Line %s did not include both units"
                    raise core.InvalidState(msg % index)
                interaction = data[(row[0], row[1])]
                interaction['unit_id_1'] = row[0]
                interaction['unit_id_2'] = row[1]
                interaction['f_crossing'] = int(row[3])
                interaction['pdb_id'] = pdb
                # added 2024-07-19 to distinguish old matlab annotations from new python ones
                interaction['program'] = 'matlab'

                family = row[2].strip()
                inter_type = self.interaction_type(family)
                if inter_type:
                    interaction[inter_type] = family

        key = op.itemgetter('unit_id_1', 'unit_id_2')
        return sorted(data.values(), key=key)

    def data(self, pdb, **kwargs):
        """Compute the interaction annotations for a pdb file.

        :pdb: The pdb id to process.
        :kwargs: Keyword arguments.
        :returns: The interaction annotations.
        """
        mlab = matlab.Matlab(str(self.config['locations']['fr3d_root']))

        self.logger.info('Running matlab on %s', pdb)
        ifn, status, err_msg = mlab.loadInteractions(pdb, nout=3)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn, pdb)
            os.remove(ifn)
            if len(data) == 0:
                data = {}
                data['pdb_id'] = pdb
                data['unit_id_1'] = 'placeholder'
                data['unit_id_2'] = 'placeholder'
                self.logger.info('PDB file %s has no interactions, save a placeholder' % pdb)
                return [data]
            return data
        elif status == 2:
            data = {}
            data['pdb_id'] = pdb                          # key is column name, value goes in the column
            data['unit_id_1'] = 'placeholder'
            data['unit_id_2'] = 'placeholder'
            self.logger.info('PDB file %s has no nucleotides, save a placeholder' % pdb)
            return [data]
        raise core.InvalidState('Matlab error code %i when analyzing %s' % (status, pdb))
