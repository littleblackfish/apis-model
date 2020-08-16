#!/bin/bash

DATADIR="../data"
MODELDIR="../models"

SPECIES="Apis_mellifera"

#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/exon/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/exon/ -b 1024 -m models.json -r 500
#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/exon/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/exon/ -b 1024 -m models.json -r 300

#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/genome/r300_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_two.json -r 300
#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_two.json -r 500

#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 128 -m models_two.json -r 300 -t 0.75 
#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 128 -m models_two.json -r 500 -t 0.75 

#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome_complete/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome_complete/ -b 128 -m models_two.json -r 300 -t 0.75 
#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome_complete/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome_complete/ -b 128 -m models_two.json -r 500 -t 0.75 

#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 512 -m models_one.json -r 500
#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_one.json -r 500

#../scripts/train_convo_cv_negative.py -d ${DATADIR}/${SPECIES}/exome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/exome/ -b 1024 -m models_one.json -r 500

SPECIES="Bombus_terrestris"

#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_one.json -r 500

#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 512 -m models_one.json -r 500 -t 0.75 

#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 512 -m models_one.json -r 500

#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_one.json -r 500

#../scripts/train_convo_cv_negative.py -d ${DATADIR}/${SPECIES}/exome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/exome/ -b 1024 -m models_one.json -r 500

SPECIES="Nasonia_vitripennis"

#../scripts/train_convo_cv.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_one.json -r 500

#../scripts/train_convo_cv_methylome.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 512 -m models_one.json -r 500 -t 0.75 

#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/methylome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/methylome/ -b 512 -m models_one.json -r 500

#../scripts/train_convo_cv_full.py -d ${DATADIR}/${SPECIES}/genome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/genome/ -b 1024 -m models_one.json -r 500

#../scripts/train_convo_cv_negative.py -d ${DATADIR}/${SPECIES}/exome/r500_onehot.cdf  -o ${MODELDIR}/${SPECIES}/exome/ -b 1024 -m models_one.json -r 500
