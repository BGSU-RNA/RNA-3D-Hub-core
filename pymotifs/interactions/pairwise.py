"""
Classes for loading pairwise interactions produced by FR3D into the RNA 3D Hub
database.
"""

import re
import os
import csv
import collections as coll

import pymotifs.core as core

from pymotifs.mat_files import Loader as MatLoader
from pymotifs.units.info import Loader as InfoLoader
from pymotifs.pdbs.infor import Loader as PdbLoader

from pymotifs.models import UnitPairsInteractions as Interaction

IGNORE = ['perp', 'nbif', 'bif']


class Loader(core.SimpleLoader):
    """A loader to generate and import the interaction annotations for
    structures.
    """

    allow_no_data = True

    dependencies = set([MatLoader, InfoLoader, PdbLoader])

    def query(self, session, pdb):
        """Create a query to access interaction data for the given pdb.

        :session: The database session to use.
        :pdb: The pdb id
        :returns: A query to get interaction data.
        """
        return session.query(Interaction).filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        """Compute the interaction annotations for a pdb file.

        :pdb: The pdb id to process.
        :kwargs: Keyword arguments.
        :returns: The interaction annotations.
        """
        matlab = core.Matlab(str(self.config['locations']['fr3d_root']))

        self.logger.info('Running matlab on %s', pdb)
        ifn, status, err_msg = matlab.loadInteractions(pdb, nout=3)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn, pdb)
            os.remove(ifn)
            return [Interaction(**d) for d in data]
        elif status == 2:
            raise core.Skip('Pdb file %s has no nucleotides' % pdb)
        raise core.InvalidState('Matlab error code %i when analyzing %s' %
                                status, pdb)

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
        else:
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
            for row in csv.reader(raw, delimiter=',', quotechar='"'):
                interaction = data[(row[0], row[1])]
                interaction['unit1_id'] = row[0]
                interaction['unit2_id'] = row[1]
                interaction['f_crossing'] = int(row[3])
                interaction['pdb_id'] = pdb

                family = row[2].strip()
                inter_type = self.interaction_type(family)
                if inter_type:
                    interaction[inter_type] = family

        return data.values()
