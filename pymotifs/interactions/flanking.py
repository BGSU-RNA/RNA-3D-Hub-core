"""
Import all flanking interactions.
Runs fr3d on the given files to determine all flankng interactions and then
imports them into the database.
"""

import os
import csv


from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod

from pymotifs.mat_files import Loader as MatLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader


class Loader(core.SimpleLoader):
    """
    A loader to generate and import the interaction annotations for
    structures.
    """

    allow_no_data = True

    dependencies = set([MatLoader, UnitLoader, PdbLoader])

    @property
    def table(self):
        return mod.UnitPairsFlanking


    def query(self, session, pdb):
        """
        Create a query to access interaction data for the given pdb.
        This is how it knows which pdb files already have data.
        However, pdb files without any flanking interactions have no data,
        so they get read week after week, wasting over 1 hour as of 2023-12-06.

        :session: The database session to use.
        :pdb: The pdb id
        :returns: A query to get interaction data.
        """
        return session.query(mod.UnitPairsFlanking).filter_by(pdb_id=pdb)


    def to_process(self, pdbs, **kwargs):
        """
        Query to find pdb ids with 'placeholder' unit ids, which means that the
        structure was processed to find pairwise interactions, none were found,
        and now we can avoid checking again.
        """

        # prevent trying to create Matlab annotations when filling in DNA releases
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        self.logger.info("nr_molecule_parent_current: %s" % nr_molecule_parent_current)

        if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
            raise core.Skip("interactions.flanking does not annotate DNA structures")

        # find all unique pdb ids that have been processed and have a flanking interaction or have a placeholder
        with self.session() as session:
            query = session.query(mod.UnitPairsFlanking.pdb_id).\
                distinct()
            done = set()
            for result in query:
                done.add(result.pdb_id)

        self.logger.info('Found %d pdbs with matlab flanking interactions' % len(done))

        needed = sorted(list(set(pdbs)-done))

        self.logger.info('Found %d pdbs that need to be processed for flanking interactions' % len(needed))

        if len(needed) == 0:
            raise core.Skip("All pdbs have been processed for flanking interactions")

        # remove the pdb ids just found from the list of those that need to be processed
        return needed


    def parse(self, filename, pdb):
        """
        Reads a csv file that lists pairs of unit ids making the flanking interaction,
        imports all interactions, deletes the file when
        done to avoid stale data and free up disk space

        :filename: The input filename.
        :pdb: The pdb id.
        :returns: A list of Interaction objects.
        """

        data = []
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for index, row in enumerate(reader):
                if not row[0] or not row[1]:
                    msg = "Line %s of %s did not include both units"
                    raise core.InvalidState(msg % (index,filename))
                data.append({'unit_id_1':row[0], 'unit_id_2':row[1], 'flanking' : int(row[2]), 'pdb_id' : pdb})

        return data


    def data(self, pdb, **kwargs):
        """
        Load and then write the the database, all flanking interactions for a pdb file.

        When the PDB file has no nucleotides or no flanking interactions, write a placeholder
        so it will not be processed again.

        :pdb: The pdb id to process.
        :kwargs: Keyword arguments.
        :returns: The interaction annotations.
        """

        # temporary code to add placeholders to all files with no flanking data
        # 2023-12-06
        # When DNA structures are added and make flanking interactions with RNA,
        # remove all placeholders, check each file once more, then run this code again.
        # data = []
        # data.append({'unit_id_1':'placeholder', 'unit_id_2':'placeholder', 'flanking' : 0, 'pdb_id' : pdb})
        # return data

        mlab = matlab.Matlab(str(self.config['locations']['fr3d_root']))  # connect to Matlab
                                                                                #### https://github.com/BGSU-RNA/RNA-3D-Hub-core/blob/09d1044cd30bd396e701d6eb91b8eef75e78b1d4/conf/bootstrap.json.txt#L31
        self.logger.info('Running matlab on %s', pdb)
        ifn, status, err_msg = mlab.loadFlankings(pdb, nout=3)            # Matlab loads .mat file for this pdb and returns flanking pair list
        # print(ifn)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn, pdb)
            os.remove(ifn)
            return data
        elif status == 2:
            self.logger.info('%s has no nucleotides, adding placeholder value', pdb)
            data = {'unit_id_1':'placeholder', 'unit_id_2':'placeholder', 'flanking' : 0, 'pdb_id' : pdb}
            return data
            #raise core.Skip('PDB file %s has no nucleotides' % pdb)
        elif status == 3:
            self.logger.info('%s has no flanking interactions, adding placeholder value', pdb)
            data = {'unit_id_1':'placeholder', 'unit_id_2':'placeholder', 'flanking' : 0, 'pdb_id' : pdb}
            return data
            # raise core.Skip('PDB file %s has no flanking interactions' % pdb)
        raise core.InvalidState('Matlab error code %i when analyzing %s' % status, pdb)
