#!/bin/bash

# Hard mask a soft masked fasta file

sed 's/[actg]/N/g' $1.fa > ${1}_masked.fa
