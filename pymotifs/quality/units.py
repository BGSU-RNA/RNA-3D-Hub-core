"""Parse and store validation reports for unit level information. This will
process the downloaded validation reports and use them to populate the
units_quality table.
"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.quality import utils as qual

from pymotifs.units.info import Loader as InfoLoader
from pymotifs.quality.download import Loader as Downloader


class Loader(core.SimpleLoader):
    """The loader to fetch and store quality data for structures.
    """
    dependencies = set([InfoLoader, Downloader])

    allow_no_data = True  # don't recompute just because there is no data
    mark = True           # note each pdb when it is processed

    def to_process(self, pdbs, **kwargs):
        """Compute the PDBs to process. These are only the PDB's that have
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
            raise core.Skip("No PDBs to process for quality.units")

        return sorted(pdbs_to_process)


    def filename(self, pdb):
        """Get the filename where the validation report of the PDB is stored.

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

    def query(self, session, pdb):
        """Generate a query to find all entries in units_quality for the given
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
        return session.query(mod.UnitQuality).\
            join(mod.UnitInfo,
                 mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
            filter(mod.UnitInfo.pdb_id == pdb)

    def as_quality(self, entry):
        """Convert an entry from the parser into a form suitable for writing to
        the units_quality table. Since some entries from the parser expand to
        more than one unit due to symmetry operators this will produce an
        iterator that may have more than 1 value.

        Parameters
        ----------
        entry : dict
            A dictionary from `Parser.nts`.

        Returns
        -------
        entry : dict
            A dictionary of 'unit_id', 'real_space_r', 'density_correlation',
            'real_space_r_z_score'.
        """
        return mod.UnitQuality(
            unit_id=entry['id'],
            real_space_r=entry.get('real_space_r', None),
            real_space_r_z_score=entry.get('real_space_r_z_score', None),
            density_correlation=entry.get('density_correlation', None),
            rscc=entry.get('rscc', None),
            rna_score=entry.get('rna_score', None),
        )

    def parse(self, filename, mapping):
        """Actually the parse the file and map quality data to unit ids.

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
            return map(self.as_quality, parser.nts(mapping))

    def data(self, pdb, **kwargs):
        """Compute the quality assignments for residues in the structure. This
        will fetch the validation report from PDB and convert the entries there
        into forms suitable to write to the database. If the report has no RSR
        or DCC data then a `core.Skip` exception will be raised.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
        data : iterable
            An iterable of a quality assignments to store in the database.
        """
        util = qual.Utils(self.config, self.session)
        mapping = util.unit_mapping(pdb)
        return self.parse(self.filename(pdb), mapping)
