"""
This is a stage to load all clash data into our database. This does not use the
atom style unit ids because we don't have any entries about this so far, but it
is simple enough to create them.
Database table is unit_clashes
"""

import datetime
import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.quality import utils as qual

from pymotifs.units.info import Loader as InfoLoader
from pymotifs.quality.download import Loader as Downloader

from fr3d.unit_ids import decode


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader, Downloader])

    allow_no_data = True  # don't recompute just because there is no data
    mark = True           # note each pdb when it is processed

    def to_process(self, pdbs, **kwargs):
        """
        Compute the PDBs to process. These are only the PDB's that have
        stored validation files and have validation data.

        Parameters
        ----------
        pdbs : list
            List of PDB ids to consider.

        Returns
        -------
        pdbs : list
            List of PDB's that have validation data to process.
        """

        # search filenames to see what pdbs have quality data
        # known = set(self._create(qual.Utils).known(has_data=True))

        # query pdb_analysis_status to see what pdbs have been processed
        with self.session() as session:
            query = session.query(mod.PdbAnalysisStatus.pdb_id).\
                filter(mod.PdbAnalysisStatus.stage == 'quality.clashes').\
                distinct()

            pdbs_processed = set()
            for result in query:
                pdbs_processed.add(result.pdb_id)

        pdbs_to_process = set(pdbs) - pdbs_processed

        if len(pdbs_to_process) == 0:
            raise core.Skip("No PDBs to process for clashes")

        return sorted(pdbs_to_process)

    def query(self, session, pdb):
        """
        Generate a query to find all entries in unit_clashes for the given
        PDB id.

        Attributes
        ----------
        session : Session
            The `Session` to use.

        pdb : str
            The PDB id to use.

        Returns
        -------
        query : Query
            Returns an SqlAlchemy query for all entires in units_quality with
            the given PDB id.
        """
        return session.query(mod.UnitClashes).\
            join(mod.UnitInfo,
                 mod.UnitInfo.unit_id == mod.UnitClashes.unit_id_1).\
            filter(mod.UnitInfo.pdb_id == pdb)

    def filename(self, pdb):
        """
        Get the filename where the validation report of the PDB is stored.

        Parameters
        ----------
        pdb : str
            The PDB id string

        Returns
        -------
        filename : str
            The path to validation report for this file.
        """
        return self._create(qual.Utils).filename(pdb)


    def compatible_units(self, unit_id1, unit_id2):
        unit1 = decode(unit_id1)
        unit2 = decode(unit_id2)
        return unit1['model'] == unit2['model']

        # Omit the line about symmetry; units with different
        # symmetries can be in the same structure and can clash
        # return unit1['model'] == unit2['model'] and \
        #     unit1['symmetry'] == unit2['symmetry']


    def as_clash(self, data):
        for unit_id1, unit_id2 in zip(*data['unit_ids']):
            if not self.compatible_units(unit_id1, unit_id2):
                raise core.InvalidState("Incompatible unit pair %s, %s" %
                                        (unit_id1, unit_id2))

            yield mod.UnitClashes(
                unit_id_1=unit_id1,
                unit_id_2=unit_id2,
                atom_name_1=data['atoms'][0],
                atom_name_2=data['atoms'][1],
                distance=data['distance'],
                magnitude=data['magnitude'],
            )

    def parse(self, filename, mapping):
        """
        Actually parse the file and map quality data to unit ids.

        Parameters
        ----------
        filename : str
            The filename to parse. The filename
        mapping : dict
            Mapping from unit keys to unit ids.

        Raises
        ------
        core.Skip
            Raised if there is not RSR and no DCC data.

        Returns
        -------
        An iterable of all unit level quality data.
        """
        with open(filename, 'rb') as raw:
            parser = qual.Parser(raw.read())
            try:
                for data in parser.clashes(mapping):
                    for clash in self.as_clash(data):
                        yield clash
            except Exception as err:
                self.logger.exception(err)
                raise core.Skip("Could not load clashes")

    def data(self, pdb, **kwargs):

        util = qual.Utils(self.config, self.session)

        # get a dictionary from simplified 5-element keys to unit ids
        # it's basically chain,residue number getting mapped to unit id
        # but it allows for model to be specified, or maybe alt id, insertion code
        mapping = util.unit_mapping(pdb)

        data = self.parse(self.filename(pdb), mapping)

        # data is an iterable of unit_clashes objects
        # those have two unit ids, atoms, distance, and magnitude

        if not data:
            raise core.Skip("No clash data found for %s" % pdb)

        return data
