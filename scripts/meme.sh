#!/usr/bin/env bash

cat positive.fa negative.fa > all.fa
fasta-get-markov -dna -m 6 all.fa all.markov
for s in {0..99}
do
  meme -dna -revcomp -objfun se -mod zoops \
  -nmotifs 10 -evt 1e-4 \
  -bfile all.markov \
  -maxw 20 \
  -seed $s \
  -neg negative.fa \
  -oc $s \
  -p 8 \
  positive.fa
done

#meme2meme -bg all.ma^Cmeme.txt > meme.out
#fimo -bfile all.markov --qv-thresh --thresh 1e-2 --max-stored-scores 10000000 meme.out all.fa
