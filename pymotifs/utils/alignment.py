import os
import logging
import shutil
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline as Clustal


logger = logging.getLogger(__name__)


def align(data):
    tmpdir = tempfile.mkdtemp()
    infile = os.path.join(tmpdir, "input.fasta")
    outfile = os.path.join(tmpdir, "output.aln")

    sequences = []
    for index, entry in enumerate(data):
        id = "seq-%s" % index
        sequences.append(SeqRecord(Seq(entry['sequence']), id=id))

    SeqIO.write(sequences, infile, "fasta")

    aligner = Clustal('clustalw2', INFILE=infile, OUTFILE=outfile)
    aligner()

    alignment = AlignIO.read(outfile, "clustal")
    shutil.rmtree(tmpdir)

    mapping = []
    indexes = [0] * len(data)
    columns = len(alignment[0])
    for col in range(columns):
        column = alignment[:, col]
        ids = []
        for index, entry in enumerate(data):
            seq_id = None
            if column[index] != '-':
                seq_id = entry['ids'][indexes[index]]
                indexes[index] += 1
            ids.append(seq_id)

        mapping.append(ids)

    return mapping
