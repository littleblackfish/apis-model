#!/usr/bin/env python
# coding: utf-8

from utils import *

import pandas as pd

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate, LeaveOneGroupOut
from sklearn.metrics import roc_curve, precision_recall_curve, auc

from argparse import ArgumentParser
from os.path import basename, splitext
import pickle

def pr_auc(model, X_test, y_test) :
    y_pred = model.predict_proba(X_test)[:,1]
    precision, recall, threshold = precision_recall_curve(y_test, y_pred)
    return auc(recall, precision)

parser = ArgumentParser()
parser.add_argument("-d", dest="data", required=True)
parser.add_argument("-c", dest='c', type=float, default=0.01)
parser.add_argument("-o", dest="out_file", default=None)
parser.add_argument("-n", dest="njobs", type=int, default=2)

args = parser.parse_args()

# Load the data

print (f'Reading {args.data}')
data = pd.read_pickle(args.data)

print (f'{len(data)} samples, {data.shape[1]} features.')
print (f'{data.attrs["n_kmers"]} kmer features.')

data = data[ (data.methylatedin == 0) | (data.methylated_ratio >0.75)]

seqid = data.index.get_level_values('seqid')

y = (data.methylatedin>1) & (data.methylated_ratio>0.75)

print (f'{len(data)} samples after filtering.')
print (f'{y.sum()} positive samples ({100*y.sum()/len(y):.1f}%).')

X = data.iloc[:,:data.attrs['n_kmers']]


model = make_pipeline(StandardScaler(), 
                      LogisticRegression(C=args.c,
                                         class_weight='balanced',
                                         penalty='l2', 
                                         solver='liblinear',
                                         verbose=0,
                                         max_iter = 1000))

params = model['logisticregression'].get_params()
print(f"Fitting {params['penalty']} penalized logistic regression with C={params['C']}.")

scores_l2 = cross_validate(model, X, y, 
                        cv=LeaveOneGroupOut().split(X,y,seqid),
                        scoring=pr_auc, 
                        n_jobs=args.njobs,
                        return_estimator=True,
                        pre_dispatch = 'n_jobs',
                        verbose = 1)

mean_score = scores_l2['test_score'].mean()  
std_score = scores_l2['test_score'].std()  

print(f'Mean score: {mean_score:.2f} (std {std_score:.2f})')

model = make_pipeline(StandardScaler(), 
                      LogisticRegression(C=args.c,
                                         class_weight='balanced',
                                         penalty='l1', 
                                         solver='liblinear',
                                         verbose=0,
                                         max_iter = 1000))

params = model['logisticregression'].get_params()
print(f"Fitting {params['penalty']} penalized logistic regression with C={params['C']}.")

scores_l1 = cross_validate(model, X, y, 
                        cv=LeaveOneGroupOut().split(X,y,seqid),
                        scoring=pr_auc, 
                        n_jobs=args.njobs,
                        return_estimator=True,
                        pre_dispatch = 'n_jobs',
                        verbose = 1)

mean_score = scores_l1['test_score'].mean()  
std_score = scores_l1['test_score'].std()  

print(f'Mean score: {mean_score:.2f} (std {std_score:.2f})')

scores = dict( features=X.columns.to_list(), l1=scores_l1, l2=scores_l2 )

if args.out_file is None : 
    root, ext = splitext(basename(args.data))
    outf = f'{root}_logistic.pkl'
else:
    outf = args.out_file

print(f'Writing model to {outf}')

with open(outf, 'wb') as f :
    pickle.dump(scores, f)
