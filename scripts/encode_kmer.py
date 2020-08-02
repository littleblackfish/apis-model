#!/usr/bin/env python
# coding: utf-8

from utils import *
import json
import xarray as xr
import pandas as pd
from argparse import ArgumentParser, FileType
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import numpy as np
from os.path import splitext,basename

parser = ArgumentParser(description="Calculates kmer counts for sequences.")
parser.add_argument( '-i',    dest="fasta_files", nargs='*', type=FileType('r'), help="Fasta file containing sequences", required=True)
parser.add_argument( '-o',    dest="out_file", type=FileType('wb'), help="Output file name", required=True )
parser.add_argument( '-k',    dest="k", type=int, required=True)
parser.add_argument( '-r',    dest="radius", type=int, default=0, help='Radius around center (default is whole window). Trimming windows will drop overlap features.')
parser.add_argument( '-norc', dest="reverse_complement", action="store_false", help='Do not include reverse complements')
args = parser.parse_args()

k = args.k 
mapping = {kmer:i for i, kmer in enumerate(all_kmers(k))}

sequences = list()
positions = list()
datas = list()

for file in args.fasta_files :
    print(f'Reading {file.name}')
    sequences += list(SeqIO.parse(file, "fasta"))

len_sequence = len(sequences[0])
assert len_sequence %2 ==0
radius = (len_sequence-2)//2

# Trimming arithmetic
assert args.radius <= radius
if args.radius and args.radius != radius : 
    trim = radius-args.radius
    assert len(sequences[0][trim:-trim]) == args.radius*2+2
    assert sequences[0].seq[trim:-trim][args.radius:args.radius+2].upper() == 'CG'
    radius = args.radius
else:
    trim =0

kmer_counts = np.zeros( (len(sequences), 4 ** k) , dtype="uint16")

print(f'Counting {k}mers in a radius of {radius}.')

for i, record in enumerate(tqdm(sequences)):

    position, data = parse_fasta_description(record)

    positions.append(position)
    datas.append(data)

    # Trim to radius
    sequence = record.seq.upper()
    assert len(sequence) == len_sequence
    if trim >0 : 
        sequence = sequence[trim:-trim]
    
    # Make sure there is a CpG in the center
    assert sequence[radius:radius+2] == 'CG'

    seq = str(sequence)


    for j in range(len(seq) - k + 1):
        kmer = seq[j : j + k]
        kmer_counts[i, mapping[ kmer ] ] +=1

positions = pd.MultiIndex.from_tuples(positions, names=["seqid", "start", "end"])

# Adjust positions for new raidus        
if trim>0 :
    tmp = positions.to_frame()
    tmp.start +=trim
    tmp.end -=trim 
    positions = pd.MultiIndex.from_frame(tmp)

datas = pd.DataFrame(datas, index=positions)

# Drop overlap features if we trimmed the windows
if trim>0 :
    try: 
        datas.drop(['gene_overlap', 'exon_overlap'], axis=1, inplace=True)
    except:
        # Maybe they weren't there to begin with
        pass

kmers= pd.DataFrame(kmer_counts, index=positions, columns=all_kmers(k))

assert (kmers.sum(axis=1) == (radius*2+2)-k+1).any()

if args.reverse_complement : 
    print ('Collapsing reverse complements..')

    reverse_complement = { kmer:str(Seq(kmer).reverse_complement()) for kmer in all_kmers(k) }

    keep, drop = set(), set()

    for kmer in all_kmers(k) : 
        if kmer not in drop : 
            keep.add(kmer)
            drop.add(reverse_complement[kmer])
    
    # Account for symmetric ones
    keep = keep.difference(drop)
        
    for kmer in keep :
        kmers[kmer] += kmers.pop(reverse_complement[kmer]) 


df = kmers.join(datas.astype('float32'))

print ('Sorting index..')
df.sort_index(inplace=True)

df.attrs=dict(k=args.k, r=radius, n_kmers=kmers.shape[1], input=[f.name for f in args.fasta_files])


print(f'Writing {args.out_file.name}')
df.to_pickle(args.out_file)