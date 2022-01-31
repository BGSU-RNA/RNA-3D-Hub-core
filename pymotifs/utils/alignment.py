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

def align_dna(data):
    tmpdir_dna = tempfile.mkdtemp()
    infile_dna = os.path.join(tmpdir_dna, "input_dna.fasta")
    outfile_dna = os.path.join(tmpdir_dna, "output_dna.aln")

    dna_sequences = []
    for index, entry in enumerate(data):
        id = "seq-%s" % index                                               ## EX: id = "seq-%s" % 5 ## print(id) => seq-5
        dna_sequences.append(SeqRecord(Seq(entry['sequence']), id=id))          ##

    SeqIO.write(dna_sequences, infile_dna, "fasta")

    aligner = Clustal('clustalw2', INFILE=infile_dna, OUTFILE=outfile)
    aligner()

    alignment = AlignIO.read(outfile, "clustal")
    shutil.rmtree(tmpdir_dna)

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

    def align_dna2(data):
        mapping = []
        if data[0] == data[1] & length(data[0])<20:
            for index, base in enumerate(data[0]):
                mapping.append(index)
    return mapping



