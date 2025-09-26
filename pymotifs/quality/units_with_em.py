"""
Parse and store validation reports for unit level information.
This will process the downloaded validation reports and use them to populate the
units_quality table.

It is not clear how clash_sum and clash_count get written to the unit_quality table;
that is not done by this stage.
"""

import os

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.quality import utils as qual  # location of the parser for xml files

from pymotifs.units.info import Loader as InfoLoader
from pymotifs.quality.download import Loader as Downloader


class Loader(core.SimpleLoader):
    """
    The loader to fetch and store quality data for structures.
    """
    dependencies = set([InfoLoader, Downloader])

    allow_no_data = True  # don't recompute just because there is no data
    mark = True           # note each pdb when it is processed

    use_marks = False     # do not use the pdb_analysis_status table entry
    merge_data = True     # allow overwriting rows with new/additional data

    testing = False       # print information, do not write to database
    fixing = False         # special code, see below

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

        if len(pdbs) == 1:
            return pdbs

        # search filenames to see what pdbs have quality data
        # known = set(self._create(qual.Utils).known(has_data=True))

        # query pdb_analysis_status to see what pdbs have been processed
        # with self.session() as session:
        #     query = session.query(mod.PdbAnalysisStatus.pdb_id).\
        #         filter(mod.PdbAnalysisStatus.stage == 'quality.units').\
        #         distinct()
        #     pdbs_processed = set([row.pdb_id for row in query])

        if self.fixing:
            # query to find EM structures, write data only for EM, Q_score, residue_inclusion
            with self.session() as session:
                query = session.query(mod.PdbInfo.pdb_id).\
                        filter(mod.PdbInfo.experimental_technique == "ELECTRON MICROSCOPY").\
                        distinct()
                em_pdbs = set([row.pdb_id for row in query])

            self.logger.info("Found %d em pdb_ids to process" % len(em_pdbs))
            print("Found %d em pdb_ids to process" % len(em_pdbs))

            pdbs = set(pdbs) & em_pdbs

            # query unit_quality table to find unit_ids we can skip
            with self.session() as session:
                query = session.query(mod.UnitInfo.pdb_id).\
                        join(mod.UnitQuality, mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
                        filter(mod.UnitQuality.Q_score.isnot(None)).\
                        filter(mod.UnitInfo.pdb_id.in_(pdbs)).\
                        distinct()
                pdbs_with_q_score = set([row.pdb_id for row in query])

            # focus on em structures that need Q_score data added
            pdbs_to_process = set(pdbs) - pdbs_with_q_score

            self.logger.info("Found %d em pdb_ids without any Q_score" % len(pdbs_to_process))
            print("Found %d em pdb_ids without any Q_score" % len(pdbs_to_process))

        else:
            # production
            # query unit_quality table to see what pdbs have data there
            with self.session() as session:
                query = session.query(mod.UnitInfo.pdb_id).\
                        join(mod.UnitQuality, mod.UnitQuality.unit_id == mod.UnitInfo.unit_id).\
                        distinct()
                pdbs_with_data = set([row.pdb_id for row in query])
            pdbs_to_process = set(pdbs) - pdbs_with_data

        # print("Found %d pdb_ids to process" % len(pdbs_to_process))

        if len(pdbs_to_process) == 0:
            raise core.Skip("No PDBs to process for quality.units")

        return sorted(pdbs_to_process)


    def filename(self, pdb):
        """
        Get the filename where the validation report of the PDB is stored.
        /usr/local/pipeline/hub-core/validation-reports/9ICY.xml.gz

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
        """
        Generate a query to find all entries in units_quality for the given
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

        if self.fixing:
            # return a query that wants non-NULL values in Q_score
            return session.query(mod.UnitQuality).\
                join(mod.UnitInfo,
                    mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
                filter(mod.UnitQuality.Q_score.isnot(None)).\
                filter(mod.UnitInfo.pdb_id == pdb)
        else:
            return session.query(mod.UnitQuality).\
                join(mod.UnitInfo,
                    mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
                filter(mod.UnitInfo.pdb_id == pdb)

    def as_quality(self, entry):
        """
        Convert an entry from the parser into a form suitable for writing to
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

        # self.logger.info('%s Q_score %10.4f residue_inclusion %10.4f' % (entry['id'],entry.get('Q_score', 0.0),entry.get('residue_inclusion', 0.0)))

        if self.testing:
            print("entry is %s" % entry)
            # self.logger.info("entry is %s" % entry)
        elif self.fixing:
            # when filling in data for cryo-em structures, only return em data
            if entry.get('Q_score', None):
                print('%s has Q_score %s and residue_inclusion %s' % (entry['id'],entry.get('Q_score', None),entry.get('residue_inclusion', None)))
            return mod.UnitQuality(
                unit_id=entry['id'],
                Q_score=entry.get('Q_score', None),
                residue_inclusion=entry.get('residue_inclusion', None)
            )
        else:
            # production
            return mod.UnitQuality(
                unit_id=entry['id'],
                real_space_r=entry.get('real_space_r', None),
                real_space_r_z_score=entry.get('real_space_r_z_score', None),
                density_correlation=entry.get('density_correlation', None),
                rscc=entry.get('rscc', None),
                rna_score=entry.get('rna_score', None),
                Q_score=entry.get('Q_score', None),
                residue_inclusion=entry.get('residue_inclusion', None)
            )

    def parse(self, filename, mapping):
        """
        Parse the file and map quality data to unit ids.

        Parameters
        ----------
        filename : str
            The filename to parse. The filename
        mapping : dict
            Mapping from unit keys to unit ids.

        Returns
        -------
        An iterable of all unit level quality data.
        """
        with open(filename, 'rb') as raw:
            parser = qual.Parser(raw.read())
            return map(self.as_quality, parser.nts(mapping))

    def data(self, pdb, **kwargs):
        """
        Compute the quality assignments for residues in the structure.
        This will read the validation report from PDB and convert the entries there
        into forms suitable to write to the database.
        If the report has no RSR
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

        filename = self.filename(pdb)
        if not os.path.exists(filename):
            raise core.Skip("No quality file downloaded for %s" % pdb)

        util = qual.Utils(self.config, self.session)
        mapping = util.unit_mapping(pdb)

        data = list(self.parse(self.filename(pdb), mapping))

        self.logger.info('Found %4d unit ids to add to unit_quality' % len(data))

        return data
