#!/bin/bash

#Convert soft masked fasta file to all uppercase

sed 's/[actg]/\U&/g' $1.fa >$1_upper.fa
