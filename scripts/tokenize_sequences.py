#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from utils import *
from tqdm import tqdm

if __name__ == "__main__":

    parser = ArgumentParser(
        description="Slides a window of size k with given stride and offset over the sequences in a fasta file. Output good for GloVe embedding."
    )
    parser.add_argument(
        "-f", dest="fasta_file", help="Fasta file containing sequences", required=True
    )
    parser.add_argument("-k", dest="k", help="Word size", type=int, required=True)
    parser.add_argument("--stride", dest="stride", help="Stride", type=int, default=1)
    parser.add_argument("--offset", dest="offset", help="Offset", type=int, default=0)
    parser.add_argument(
        "-rc",
        dest="reverse_complement",
        help="Also co reverse complements",
        action="store_true",
    )
    args = parser.parse_args()

    # GloVe takes each document in a single line.
    # Spaces separate words, newlines separate documents
    # In our case, a word is a kmer, a document is a sequence.

    for record in tqdm(SeqIO.parse(args.fasta_file, "fasta")):
        sequence = record.seq.upper()
        for kmer in kmer_iter(
            str(sequence), k=args.k, stride=args.stride, offset=args.offset
        ):
            print(kmer, end=" ")
        print("")
        if args.reverse_complement:
            for kmer in kmer_iter(
                str(sequence.reverse_complement()),
                k=args.k,
                stride=args.stride,
                offset=args.offset,
            ):
                print(kmer, end=" ")
            print("")
