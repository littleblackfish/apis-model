#!/bin/bash

# Corpus body
CORPUS_FILE=$1
# Word size
K=$2
# Embedding dimensionality 
D=$3
# Output dir
OUTDIR=$4

mkdir -p ${OUTDIR}

# Filelist 
VOCAB_FILE=${OUTDIR}/vocab_${K}mer.txt
COOCCURRENCE_FILE=${OUTDIR}/cooccurrence_${K}mer.bin
COOCCURRENCE_SHUF_FILE=${OUTDIR}/cooccurrence_${K}mer.shuf.bin
SAVE_FILE=${OUTDIR}/vectors_${K}mer_${D}d
LOG_FILE=${OUTDIR}/glove_${K}mer_${D}d.log

VERBOSE=2
NUM_THREADS=8
MEMORY=8.0

# Number of possible kmers
MAX_VOCAB=$(( 4 ** $K ))
WINDOW_SIZE=$(( 5 * $K ))

MAX_ITER=1000

X_MAX=10

# Count kmer usage
[[ $VOCAB_FILE -nt $CORPUS_FILE ]] || vocab_count -max-vocab $MAX_VOCAB -verbose $VERBOSE < $CORPUS_FILE > $VOCAB_FILE
# Count kmer cooccurrances
[[ $COOCCURRENCE_FILE -nt $CORPUS_FILE ]] || cooccur -memory $MEMORY -vocab-file $VOCAB_FILE -verbose $VERBOSE -window-size $WINDOW_SIZE < $CORPUS_FILE > $COOCCURRENCE_FILE
# Shuffle cooccurrances
[[ $COOCCURRENCE_SHUF_FILE -nt $COOCCURRENCE_FILE ]] || shuffle -memory $MEMORY -verbose $VERBOSE < $COOCCURRENCE_FILE > $COOCCURRENCE_SHUF_FILE
# Train embedding
glove -save-file $SAVE_FILE -threads $NUM_THREADS -input-file $COOCCURRENCE_SHUF_FILE -x-max $X_MAX -iter $MAX_ITER -vector-size $D -binary 0 -vocab-file $VOCAB_FILE -verbose $VERBOSE 2>&1 | tee $LOG_FILE
