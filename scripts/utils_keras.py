#!/usr/bin/env python
# coding: utf-8

import numpy as np
from tensorflow import keras
from utils import *

## Generate a hybrid keras model with embedding, convolutional, recurrent and dense layers

def convo_model(
    input_shape,
    name="convo_model",
    conv_nfilters=0,
    conv_ksize=16,
    conv_psize=8,
    dense_units=64,
    dropout=0,
    lrate=1e-3,
    verbose=False
):

    model = keras.Sequential(name=name)

    ## Input layer, mainly to specify shape
    model.add(keras.layers.InputLayer(input_shape))
    
    if conv_nfilters>0 :
        model.add(
            keras.layers.Conv1D(
                filters=conv_nfilters,
                kernel_size=conv_ksize,
                strides=1,
                activation="relu"
            )
        )
        if conv_psize > 0 :
            model.add(keras.layers.MaxPool1D(pool_size=conv_psize))


    model.add(keras.layers.Flatten())
    
    if dense_units>0 :
        model.add(keras.layers.Dense(dense_units, activation="relu", name=f"dense"))
    
    if dropout > 0 : 
        model.add(keras.layers.Dropout(dropout, name=f"Dropout"))

    ## Output layer
    model.add(keras.layers.Dense(1, activation="sigmoid", name="Output"))

    model.compile(
        loss=keras.losses.BinaryCrossentropy(),
        optimizer=keras.optimizers.Adam(lrate),
        metrics=[keras.metrics.AUC(curve="PR")]
    )
    
    if verbose :
        model.summary()

    return model


## Set pretrained embedding
def set_embedding(model, mapping):

    embedding_shape = model.layers[0].get_weights()[0].shape

    # infer k
    k = int(np.log2(embedding_shape[0]) // 2)

    w = np.zeros(embedding_shape)

    for i, kmer in enumerate(all_kmers(k)):
        assert len(mapping[kmer]) == embedding_shape[1]
        w[i] = mapping[kmer]

    model.layers[0].weights[0].assign(w)


## Set output bias to log(npos/nneg)
def set_output_bias(model, y_train=None):
    num_positive = sum(y_train)
    num_negative = len(y_train) - sum(y_train)

    model.layers[-1].bias.assign([np.log(num_positive / num_negative)])

    # Return class weights for optional balancing
    return {0: 1, 1: num_negative / num_positive}

import xarray as xa
def get_data(filename, radius=None, threshold=0.75, mode='normal') :
    
    data = xa.load_dataset(filename)
    assert data.radial_pos.min().equals(- data.radial_pos.max())
    max_radius = float(data.radial_pos.max())
    assert radius is None or radius <= max_radius
    
    if mode=='normal' :
        mask = (data.methylatedin == 0) | (data.methylated_ratio > threshold)
        y = (data.methylatedin > 1) & (data.methylated_ratio > threshold)
    
    elif mode=='negative' : 
        mask = (data.methylatedin == 0) | (data.methylated_ratio < threshold)
        y = (data.methylatedin > 1) & (data.methylated_ratio < threshold)

    elif mode=='methylome' :
        mask = data.methylatedin > 1 
        y = data.methylated_ratio>threshold
    
    else:
        raise ValueError('Invalid mode, choose normal, negative, or methylome.') 
    
    seqid = data.seqid[mask]
    y = y[mask]
    X = data.onehot[mask]
    
    if radius is not None and radius < max_radius :
        X = X.where(np.abs(X.radial_pos) <= radius, drop=True) 

    assert X.position.equals(y.position)

    X = X.values
    y = y.values
    seqid = seqid.values
       
    filename = '/'.join(filename.split('/')[-2:])

    print (f'{filename}: size: {len(y)}, pos: {y.sum()}, baseline: {y.sum()/len(y):.2f}')
    
    return X,y,seqid

from tqdm.notebook import tqdm as tqdm_notebook
from glob import glob

def test_cv(model_path, data, verbose=False, batch_size=1024, tqdm=tqdm_notebook) :

    X,y,seqid = data
    
    scores = list()
    
    for path in tqdm(list(glob(f'{model_path}_NC*.h5')), leave=False) :
        
        test_seqid = path[len(f'{model_path}_'):-3 ]
               
        model = keras.models.load_model(path)
        
        test = seqid == test_seqid
        X_test = X[test]
        y_test = y[test]
        
        loss, auc = model.evaluate(X_test,y_test, batch_size=batch_size, verbose=verbose)
        
        scores.append(auc)
        
    scores = np.array(scores)
     
    print (f'mean: {scores.mean():.2f} std:{scores.std():.2f}')
    
    return np.array(scores)

def test_ensemble(model_path, data, verbose=False, batch_size=1024, tqdm=tqdm_notebook) :

    X,y = data

    scores = list()
    
    for path in tqdm(list(glob(f'{model_path}_NC*.h5')), leave=False) :
                        
        model = keras.models.load_model(path)
        
        loss, auc = model.evaluate(X,y, batch_size=batch_size, verbose=verbose)
        
        scores.append(auc)
        
    scores = np.array(scores)
     
    print (f'mean: {scores.mean():.2f} std:{scores.std():.2f}')
    
    return scores
