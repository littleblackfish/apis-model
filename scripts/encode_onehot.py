#!/usr/bin/env python
# coding: utf-8

from utils import *
import json
import numpy as np
import pandas as pd
import xarray as xr
from Bio import SeqIO
from tqdm import tqdm

from argparse import ArgumentParser,FileType

parser = ArgumentParser(description="Encodes sequences to one-hot representations.")
parser.add_argument( '-i',    dest="fasta_files", nargs='*', type=FileType('r'), help="Fasta file containing sequences", required=True)
parser.add_argument( '-o',    dest="out_file", help="Output file name", required=True )
parser.add_argument( '-norc', dest="reverse_complement", action="store_false", help='Do not include reverse complements')
args = parser.parse_args()

def onehot_encode(sequence):
    # I suspect this is marginally faster than a dict mapped solution for small alphabet but haven't tested rigorously. 
    encoded = np.zeros((len(sequence), 4), dtype=bool)
    for i, letter in enumerate(sequence.upper()):
        if letter == 'A':
            encoded[i, 0] = 1
        elif letter == 'T':
            encoded[i, 1] = 1
        elif letter == 'C':
            encoded[i, 2] = 1
        elif letter == 'G':
            encoded[i, 3] = 1
        elif letter != 'N' :
            raise ValueError('Weird ass character in sequence:', letter)

    return encoded

# Parse sequences from fasta

for file in args.fasta_files :
    print(f'Reading {file.name}')

    if args.reverse_complement :
        # It is best to have reverse complement consecutive to the sequence to make sure they end up in the same batch
        print("Parsing sequences and generating reverse complements..")

        sequences = list()

        for record in tqdm(SeqIO.parse(file, "fasta")):
            sequences.append(record)

            pos, data = parse_fasta_description(record)

            seq = record.seq.reverse_complement()

            id = f"{pos[0]}:{pos[2]}-{pos[1]}"

            description = json.dumps(data, separators=(",", ":"))

            sequences.append(
                SeqIO.SeqRecord(seq=seq, id=id, description=f"{id} {description}")
            )
    
    if not args.reverse_complement : 
        print("Parsing sequences..")
        sequences = list(tqdm(SeqIO.parse(file, "fasta")))

len_sequence = len(sequences[0])
radius = (len_sequence - 2) / 2

positions = list()
datas = list()

print("One-hot encoding...")

onehot = np.empty((len(sequences), len_sequence, 4), dtype=bool)

for i, record in enumerate(tqdm(sequences)):
    
    position, data = parse_fasta_description(record)
    positions.append(position)
    datas.append (data)
    
    onehot[i] = onehot_encode(record.seq)

positions = pd.MultiIndex.from_tuples(positions, names=["seqid", "start", "end"])

datas = pd.DataFrame(datas, index=positions).astype('float32')

radial_pos = np.concatenate([np.arange(-radius, 1, dtype=int), np.arange(radius+1, dtype=int) ])
letters = ["A", "T", "C", "G"]

onehot = xr.DataArray(
    onehot,
    dims=("position", "radial_pos", "letter"),
    coords=dict( position=positions, radial_pos=radial_pos, letter=letters ),
    name="onehot"
)

ds = xr.Dataset(dict(onehot=onehot))

ds = ds.assign( {feature:xr.DataArray(datas[feature], dims=("position")) for feature in datas })

ds.reset_index("position").to_netcdf(args.out_file)
