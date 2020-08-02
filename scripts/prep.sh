#!/bin/bash

R=$1
SCRIPT="../../scripts"

${SCRIPT}/extract_windows.py -in ../robust_CpG.pkl -genome ../Amel_HAv3.1.fa -r $R -r_exclude 50 -exon -chromosome -o robust_$R.fa
${SCRIPT}/extract_windows.py -in ../never_CpG.pkl -genome ../Amel_HAv3.1.fa -r $R -r_exclude 50 -exon -chromosome -o never_$R.fa
${SCRIPT}/upper.sh robust_$R
${SCRIPT}/upper.sh never_$R
${SCRIPT}/encode_sequences_kmer.py -i robust_${R}_upper.fa -o robust_kmer_${R}.cdf 
${SCRIPT}/encode_sequences_kmer.py -i never_${R}_upper.fa -o never_kmer_${R}.cdf 
${SCRIPT}/merge_pos_neg.py -p robust_kmer_${R}.cdf -n never_kmer_${R}.cdf -o kmer_${R}.cdf
