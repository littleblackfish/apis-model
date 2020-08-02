#!/usr/bin/env python
# coding: utf-8

from utils import *
import json,re
import xarray as xr
import pandas as pd
from argparse import ArgumentParser, FileType
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import numpy as np
from os.path import splitext,basename

parser = ArgumentParser(description="Calculates dreme motif counts for sequences.")
parser.add_argument( '-i',    dest="fasta_files", nargs='*', type=FileType('r'), help="Fasta file containing sequences", required=True)
parser.add_argument( '-d',    dest="dreme_file", type=FileType('r'), help='Dreme motif file', required=True)
parser.add_argument( '-o',    dest="out_file", type=FileType('wb'), help="Output file name", required=True )
parser.add_argument( '-norc', dest="reverse_complement", action="store_false", help='Do not include reverse complements')
args = parser.parse_args()


print (f'Parsing motifs from {args.dreme_file.name}')

motifs = []

# super primitive dreme output parser

for line in args.dreme_file :
    if line[:5]=='MOTIF' :
        tmp, motif, tmp2 = line.split()
        motifs.append(motif)

# replace iupac chars with regular expression equivalents

iupac =dict(#A='A', C='C', T='T', G='G', 
           W='[AT]', S='[CG]', M='[AC]', K='[GT]', R='[AG]', Y='[CT]',
           B='[CGT]', D='[AGT]', H='[ACT]', V='[ACG]', N='[ACGT]')

re_motifs=[]

for motif in motifs:
    for key in iupac :
        motif = motif.replace(key, iupac[key])
    re_motifs.append(motif)


sequences = list()
positions = list()
datas = list()

for file in args.fasta_files :
    print(f'Reading {file.name}')
    sequences += list(SeqIO.parse(file, "fasta"))

len_sequence = len(sequences[0])

print(f'Counting motifs.')

counts = np.zeros((len(sequences), len(motifs)), dtype='uint16')
index = []

for i, record in enumerate(tqdm(sequences)) :
    
    index.append(parse_fasta_description(record)[0][:2])

    position, data = parse_fasta_description(record)

    sequence = record.seq.upper()

    positions.append(position)
    datas.append(data)

    for j, motif in enumerate(re_motifs) :
        counts[i,j] = len(re.findall(motif, str(sequence)))
        if args.reverse_complement :
            counts[i,j] += len(re.findall(motif, str(sequence.reverse_complement())))

        
positions = pd.MultiIndex.from_tuples(positions, names=["seqid", "start", "end"])

datas = pd.DataFrame(datas, index=positions)

counts= pd.DataFrame(counts, index=positions, columns=motifs)

df = counts.join(datas.astype('float32'))

print ('Sorting index..')
df.sort_index(inplace=True)

df.attrs=dict(n_kmers=len(motifs), input=[f.name for f in args.fasta_files])


print(f'Writing {args.out_file.name}')
df.to_pickle(args.out_file)