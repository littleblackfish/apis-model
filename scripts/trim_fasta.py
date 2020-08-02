#!/usr/bin/env python3

from Bio import SeqIO
from sys import argv
import argparse

parser = argparse.ArgumentParser(description='Crop all the sequences in a fasta file')
parser.add_argument('filename', help='fasta filename')
parser.add_argument('-s', type=int, default=0, help='start index, 0-based, inclusive')
parser.add_argument('-e', type=int, help='end index, 0-based, inclusive')
parser.add_argument('-o', help='Output (fasta) filename, STDOUT by default')
args = parser.parse_args()

records = [record for record in SeqIO.parse(args.filename, 'fasta')]

if args.e is None : 
    for record in records :
        record.seq=record.seq[args.s:]
else :
    for record in records :
        record.seq=record.seq[args.s:args.e+1]

if args.o is None : 
    for record in records : 
        print (record.format('fasta'))
else:
    SeqIO.write(records, args.o , 'fasta')
