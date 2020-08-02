from boilerplate import *

from utils import *
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser(description="Extract exons from a genome")
parser.add_argument(
    "-o", dest="outfile", help="Output file", default="../data/exons.fa"
)
args = parser.parse_args()

exons = list()

for mrna in representative_mrnas:
    for exon in db.children(mrna, featuretype="exon"):
        id = f"{exon.seqid}:{exon.start}-{exon.end}[{exon.strand}]"
        seq = feature_sequence(genome, exon)
        assert len(seq) == len(exon)
        exons.append(SeqRecord(seq, id))

print(f"Writing exons to {args.outfile}")
SeqIO.write(exons, args.outfile, "fasta")
