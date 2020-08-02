#!/usr/bin/env python
# coding: utf-8

from utils import *
import json
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

import pyranges as pr


from argparse import ArgumentParser

# Samples cpgs with at least radius (genomic) spacing between them
# cpgdf is any DataFrame indexed by pos

parser = ArgumentParser(description="Extracts windows around CpGs")

parser.add_argument('-i', dest='infile', required=True,
    help="Input file is a pickled DataFrame indexed by genomic position." )

parser.add_argument("-o", dest="outfile", required=True, help="Output file (fasta)")

parser.add_argument('-fa', dest="genome", required=True, 
    help="Genome file in fasta format." )

parser.add_argument('-gff', dest="gff", required=True, 
    help="Annotation file in fasta format." )

parser.add_argument( "-r", dest="radius", type=int, required=True,
    help="Radius, total winsize will be 2radius+2"
)

parser.add_argument(
    "-r_exclude", dest="overlap", type=int, default = 0,
    help="Exclusion radius; no overlap is allowed within this radius."
)

parser.add_argument(
    "-exon", dest="exon", action='store_true',
     help="Sample exclusively from exons"
)

parser.add_argument(
    "-chromosome", dest="chromosome", action='store_true',
     help="Sample exclusively from complete chromosomes"
)

# Checks if two ranges overlap
def overlap(pos1, pos2, radius):
    return pos1[0] == pos2[0] and abs(pos1[1] - pos2[1]) < 2 * radius

# Returns subset with no overlapping ranges
# Assuming poslist is sorted!
def no_overlap_subset(poslist, radius):
    spaced = [poslist[0]]

    for pos in poslist[1:]:
        if not overlap(pos, spaced[-1], radius):
            spaced.append(pos)

    return spaced

# Returns window with radius around the pos
def get_window(genome, position, radius):
    seqid, pos = position
    assert genome[seqid][pos:pos+2].upper() == 'CG'
    start = pos-radius
    end  = pos+radius+2
    seq = genome[seqid][start:end]
    return SeqIO.SeqRecord(seq=seq, id=f'{seqid}:{start}-{end}')

args = parser.parse_args()

# Load master DataFrame
print(f'Reading {args.infile}')
index = pd.read_pickle(args.infile).sort_index()

print(f"{len(index)} sites. ")

# Load genome
genome = load_genome(args.genome, upper=False)

print(f"Loading annotation from {args.gff}")
# Load annotation
annotation = pr.read_gff3(args.gff)

genes = annotation[annotation.Feature=='gene'].merge(strand=False)
exons = annotation[annotation.Feature=='exon'].merge(strand=False)

tmp = index.index.to_frame(index=False)
tmp['Start'] = tmp.pos-args.radius
tmp['End'] = tmp.pos+args.radius+2
tmp['Chromosome'] = tmp.seqid
windows = pr.PyRanges(tmp)

coverage = windows.coverage(exons)
tmp = coverage.df.set_index(['seqid', 'pos']).reindex(index.index)
index['exon_overlap'] = tmp.FractionOverlaps.round(2)

coverage = windows.coverage(genes)
tmp = coverage.df.set_index(['seqid', 'pos']).reindex(index.index)
index['gene_overlap'] = tmp.FractionOverlaps.round(2)

# Filter exonic sites
if args.exon:
    index = index[index.inexon>0]
    print(f"{len(index)} exonic sites.")

# Filter chromosomal sites
if args.chromosome:
    index = index[index.nuclear]
    print(f"{len(index)} nuclear sites.")

positions = index.index.to_list()

if args.overlap>0 :
    print("Checking overlaps.")
    positions = no_overlap_subset(positions, args.overlap)
    print(f"{len(positions)} non-overlapping sites.")
    new_index = pd.MultiIndex.from_tuples(positions, names=index.index.names)                                                                                                                              
    index = index.reindex(new_index)

print("Extracting windows..")

try :
    index.methylated_ratio = index.methylated_ratio.round(2)
except:
    index['methylated_ratio']=0

try:
    index.cond_mean_methylation = index.cond_mean_methylation.round(2)
except :
    index['cond_mean_methylation'] = 0


records = [] 
for position, row in tqdm(index[['gene_overlap', 'exon_overlap', 'inexon', 'methylatedin', 'methylated_ratio', 'cond_mean_methylation']].iterrows()) : 

    record = get_window(genome, position, args.radius)

    record.description = json.dumps(row.to_dict(),separators=(",", ":"))

    records.append(record)

# filter out short length  windows
records = filter(lambda r: len(r) == args.radius * 2 + 2, records)
# filter out windows containing Ns
records = filter(lambda r: r.seq.find("N") < 0, records)
# Write sequences
n = SeqIO.write(records, args.outfile, "fasta")

print(
    f"{n} sites written after final filtering. "
)
