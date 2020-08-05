#!/usr/bin/env python
# coding: utf-8
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
#os.environ["CUDA_VISIBLE_DEVICES"]="3"

from utils import *
from utils_keras import *
import numpy as np
import xarray as xa
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import cross_validate, LeaveOneGroupOut
from sklearn.metrics import roc_curve, precision_recall_curve, auc
import json,pickle
from os import makedirs
from os.path import join

if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-d", dest="data", required=True)
    parser.add_argument("-m", dest="models", required=True)
    parser.add_argument("-o", dest="out_dir", default='models')

    parser.add_argument('-r', dest='radius', help='Radius around center CpG', type=int, required=True)

    parser.add_argument("-b", dest="batch_size", help="Batch_size", type=int, default=64)
    parser.add_argument("-e", dest="epochs", help="Number of epochs to train", type=int, default=100)
    args = parser.parse_args()

    print (f'Reading {args.data}')
    # Load the data
    data = xa.load_dataset(args.data)

    assert data.radial_pos.min().equals(- data.radial_pos.max())
    max_radius = float(data.radial_pos.max())

    assert args.radius <= data.onehot.shape[1] 

    mask = (data.methylatedin == 0) | (data.methylated_ratio >0.75)

    seqid = data.seqid[mask]

    y = (data.methylatedin>1) & (data.methylated_ratio>0.75)
    y = y[mask]
    
    X = data.onehot[mask]

    if args.radius < max_radius :
        X = X.where(np.abs(X.radial_pos)<=args.radius, drop=True) 

    assert X.position.equals(y.position)
    
    num_positive = float(y.sum())
    num_negative = len(y) - num_positive

    print (f'{len(data.onehot)} samples.')
    print (f'{len(y)} samples after filtering.')
    print (f'{num_positive} positive samples ({100*num_positive/len(y):.1f}%).')
    
    with open(args.models) as f:
        models = json.load(f)

    for model_params in models:

        keras.backend.clear_session()
        
        model_name = model_params['name']

        log_dir = join(args.out_dir, 'log_dir', f'{model_name}_{args.radius}' )
        makedirs(log_dir, exist_ok=True)

        model_filename = join(args.out_dir, f'{model_name}_{args.radius}')

        print (f'Will write model to {model_filename}')

        model = convo_model( input_shape=X[0].shape, verbose=True, **model_params)


        callbacks = [ 
            keras.callbacks.EarlyStopping( monitor="val_loss", patience=5, restore_best_weights=True, verbose=1 ),
            keras.callbacks.ReduceLROnPlateau( monitor="val_loss", factor=0.5, patience=2, verbose=1 ),
            keras.callbacks.TensorBoard(log_dir=log_dir)
            ]

        fit_params =dict(
            callbacks=callbacks,
            batch_size=args.batch_size,
            epochs=args.epochs,
            verbose=1,
        ) 
        scores = list()

        keras.backend.clear_session()


        for train, val in LeaveOneGroupOut().split(X, y, X.seqid) :
            X_train = X.values[train]
            y_train = y.values[train]

            X_val = X.values[val]
            y_val = y.values[val]

            break
        
        model.fit(X_train, y_train, validation_data = (X_val, y_val), **fit_params) 

        model.save(f'{model_filename}_full.h5')