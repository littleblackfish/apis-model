#!/usr/bin/env python
# coding: utf-8

from utils import *
import json
import xarray as xr
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO
from tqdm import tqdm
import numpy as np


def onehot_encode(sequence):

    encoded = np.zeros((len(sequence), 4), dtype=bool)

    for i, letter in enumerate(sequence):
        if letter == "A":
            encoded[i, 0] = 1
        elif letter == "T":
            encoded[i, 1] = 1
        elif letter == "C":
            encoded[i, 2] = 1
        elif letter == "G":
            encoded[i, 3] = 1

    return encoded

parser = ArgumentParser(description="Encodes sequences in one-hot and kmer representations")
parser.add_argument(
    "-f", dest="fasta_file", help="Fasta file containing sequences", required=True
)
parser.add_argument(
    "-no_onehot", dest="onehot", help="Omit one-hot encoding", action="store_false"
)
parser.add_argument(
    "-no_kmer", dest="kmer", help="Omit kmer encoding", action="store_false"
)
parser.add_argument(
    "-rc",
    dest="reverse_complement",
    help="Introduce reverse complements.",
    action="store_true",
)
parser.add_argument("-o", dest="outfile", required=True)
args = parser.parse_args()

if not args.reverse_complement : 
    print("Reading sequences..")
    sequences = list(tqdm(SeqIO.parse(args.fasta_file, "fasta")))
else:
    # It is best to have reverse complement consecutive to the sequence to make sure they end up in the same batch
    print("Reading sequences and generating reverse complements..")

    sequences = list()

    for record in tqdm(SeqIO.parse(args.fasta_file, "fasta")):
        sequences.append(record)

        pos, data = parse_fasta_description(record)

        seq = record.seq.reverse_complement()

        id = f"{pos[0]}:{pos[2]}-{pos[1]}"

        description = json.dumps(data, separators=(",", ":"))

        sequences.append(
            SeqIO.SeqRecord(seq=seq, id=id, description=f"{id} {description}")
        )

len_sequence = len(sequences[0])
radius = (len_sequence - 2) / 2

for_k = [3, 6]

positions = list()
datas = dict()

print("Parsing positions")

for i, record in enumerate(tqdm(sequences)):
    position, data = parse_fasta_description(record)
    positions.append(position)
    datas[position] = data

positions = pd.MultiIndex.from_tuples(positions, names=["seqid", "start", "end"])

datas = pd.DataFrame.from_dict(datas, orient="index").reindex(positions)
# datas.index= positions

if args.onehot:
    print("One-hot encoding...")

    onehot = np.empty((len(sequences), len(sequences[0]), 4), dtype=bool)

    for i, record in enumerate(tqdm(sequences)):
        onehot[i] = onehot_encode(record.seq)

    onehot = xr.DataArray(
        onehot,
        dims=("position", "index", "letter"),
        coords=dict(
            position=positions,
            index=np.arange(len_sequence, dtype="int16") - int(radius),
            letter=["A", "T", "C", "G"],
        ),
        name="onehot_encoded",
    )

if args.kmer:
    print("Kmer encoding and spectrum")
    kmer_encoded = list()
    for k in for_k:
        print(f"for k={k}...")

        mapping = {kmer: i for i, kmer in enumerate(all_kmers(k))}
        bins = np.arange(4 ** k + 0.5) - 0.5

        encoded = np.empty((len(sequences), len_sequence - k + 1), dtype="uint16")
        spectrum = np.empty((len(sequences), 4 ** k), dtype="uint16")

        for i, record in enumerate(tqdm(sequences)):
            seq = str(record.seq)
            # Encode kmers to integers
            for j in range(len_sequence - k + 1):
                encoded[i, j] = mapping[seq[j : j + k]]
            # Histogram of the encoded sequence is the kmer spectrum
            spectrum[i] = np.histogram(encoded[i], bins)[0]

        kmer_encoded.append(
            xr.DataArray(
                encoded,
                dims=("position", f"{k}mer_index"),
                coords=dict(position=positions),
                attrs=dict(k=k),
                name=f"{k}mer_encoded",
            )
        )

        kmer_encoded.append(
            xr.DataArray(
                spectrum,
                dims=("position", f"{k}mer_count"),
                coords={"position": positions, f"{k}mer_count": list(all_kmers(k))},
                attrs=dict(k=k),
                name=f"{k}mer_spectrum",
            )
        )

    kmer_encoded = xr.Dataset({da.name: da for da in kmer_encoded})

master = kmer_encoded.assign(
    dict(
        onehot_encoded=onehot,
    #    mean_methylation=xr.DataArray(datas.mean_met, dims="position"),
    #    std_methylation=xr.DataArray(datas.std_met, dims="position"),
    )
)

master.reset_index("position").to_netcdf(args.outfile)
