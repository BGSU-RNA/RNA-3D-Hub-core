"""

Program for exporting data from the RNA 3D Hub database into compressed files.
Used for data exchange with NDB and for providing downloadable files.

Usage:
python PdbFileExporter.py <data_type> <path>

Example PdbFileExporter.py interactions ~/interactions.csv.gz

data types: interactions
TODO: add motifs and loops.

"""

import logging
import sys
import gzip
import tempfile
import shutil

from sqlalchemy import distinct, or_


from models import session, PairwiseInteractions, PdbUnitIdCorrespondence
from MotifAtlasBaseClass import MotifAtlasBaseClass


class PdbFileExporter(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        """
        """
        MotifAtlasBaseClass.__init__(self)
        self.output = ''

    def _get_all_pdbs_with_interactions(self):
        """
            List all pdbs with stored pairwise interaction annotations.
        """
        pdb_ids = []
        for id in session.query(distinct(PairwiseInteractions.pdb_id)).all():
            pdb_ids.append(id[0])
        return pdb_ids

    def _get_interactions(self, pdb_id):
        """
            Not all interactions are sent to the output, so need to make sure
            that all of the desired interactions are not empty.
        """
        return session.query(PairwiseInteractions).\
                       filter_by(pdb_id=pdb_id).\
                       filter(or_(PairwiseInteractions.f_lwbp   != None,
                                  PairwiseInteractions.f_stacks != None,
                                  PairwiseInteractions.f_bphs   != None)).\
                       all()

    def _get_id_correspondence(self, pdb_id):
        """
            Get a dictionary with old and new ids.
        """
        unit_id = dict()
        for ids in session.query(PdbUnitIdCorrespondence).\
                           filter_by(pdb=pdb_id).\
                           all():
            unit_id[ids.old_id] = ids.unit_id
        return unit_id

    def _format_nones(self, value):
        """
            Replace nones with empty strings so that joining doesn't fail.
        """
        return value if value else ''

    def _format_interactions(self, interaction, unit_id):
        """
            Manually create a csv-formatted string.
        """
        return '"%s"\n' % '","'.join([unit_id[interaction.iPdbSig],
                                      unit_id[interaction.jPdbSig],
                                      self._format_nones(interaction.f_lwbp),
                                      self._format_nones(interaction.f_stacks),
                                      self._format_nones(interaction.f_bphs)])

    def _create_compressed_output_file(self, f_in, output_file):
        """
            Take the file handle and send its contents to a gzipped output file.
        """
        temp_output_file = output_file + 'temp'
        f_in.seek(0)
        f_out = gzip.open(temp_output_file, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        # rename file only when it's ready
        shutil.move(temp_output_file, output_file)

    def export_interactions(self, output_file, pdb_ids=None):
        """
            Pdb files for which new-style unit ids haven't been generated are
            skipped for now.
        """
        if not pdb_ids:
            pdb_ids = self._get_all_pdbs_with_interactions()

        temp = tempfile.TemporaryFile()
        headers = ['unit_id1', 'unit_id2', 'FR3D basepair (f_lwbp)',
                   'FR3D stacking (f_stacks)', 'FR3D base phosphate (f_bphs)']
        temp.write( '"' + '","'.join(headers) + '"\n')

        for pdb_id in pdb_ids:
            logging.info('Writing out interactions for %s' % pdb_id)
            interactions = self._get_interactions(pdb_id)
            unit_id = self._get_id_correspondence(pdb_id)

            if len(unit_id) == 0:
                logging.warning('Unit ids not found for %s' % pdb_id)
                continue

            for interaction in interactions:
                temp.write(self._format_interactions(interaction, unit_id))

        self._create_compressed_output_file(temp, output_file)
        temp.close()


def main(argv):
    """
    """
    d = PdbFileExporter()
    d.start_logging()

    if argv[0] == 'interactions':
        d.export_interactions(argv[1])
    else:
        print "Unrecognized option"
        sys.exit(1)

    status = 'Data files successfully created'
    d.set_email_subject(status)
    logging.info(status)
    d.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])
