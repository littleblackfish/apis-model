#!/usr/bin/env python
# coding: utf-8

from utils import *
import json
from argparse import ArgumentParser
from Bio import SeqIO
from tqdm import tqdm
from random import sample

parser = ArgumentParser(
    description="Samples random sets of sequences from a given fasta file."
)
parser.add_argument(
    "-f", dest="fasta_file", help="Fasta file containing sequences", required=True
)
parser.add_argument(
    "-nsamples", dest="nsamples", help="number of samples", type=int, default=1
)
parser.add_argument("-nsets", dest="nsets", help="number of sets", type=int, default=1)
parser.add_argument("-filter", dest="filter", help="threshold", nargs=2, required=False)
parser.add_argument("-o", dest="outfile", help="Output filename prefix", required=True)
args = parser.parse_args()

print("Reading sequences")

if args.filter is None:
    sequences = [record for record in SeqIO.parse(args.fasta_file, "fasta")]
else:
    key = args.filter[0]
    threshold = float(args.filter[1])
    sequences = [
        record
        for record in SeqIO.parse(args.fasta_file, "fasta")
        if json.loads(record.id)[key] > threshold
    ]

print(
    f"Read {len(sequences)} sequences, sampling {args.nsets} sets of {args.nsamples}.."
)

if args.nsamples < 1:
    SeqIO.write(sequences, f"{args.outfile}.fa", "fasta")
else:
    for set in tqdm(range(args.nsets)):
        select_sequences = sample(sequences, args.nsamples)

        SeqIO.write(select_sequences, f"{args.outfile}-{set}.fa", "fasta")
