"""Utilities for training and running the model."""

#from itertools import product
import os

import numpy
import pandas as pd
import pickle
from sklearn.linear_model import RidgeCV
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler


def TransposeRnaSeqTable(rnaseq):
    """Convert the RnaSeq table, indexed by gene, to be indexed by specimen."""
    rnaseq.index = rnaseq.Gene
    return rnaseq[rnaseq.columns[2:]].T


def NormSpecimens(specimens):
    #  normed_specimens = specimens.apply(
    #      lambda specimen : specimen / numpy.linalg.norm(specimen), axis=1)
    normed_specimens = specimens.div(numpy.linalg.norm(specimens, axis=1), axis=0)
    return normed_specimens


def ApplySelectedGenes(model_dir, data):
    with open(os.path.join(model_dir, 'selected_genes.pkl'), 'rb') as f:
        selected_genes = pickle.load(f)
    data_selected = data[selected_genes]
    return data_selected


def ApplyStandardScaler(model_dir, data):
    with open(os.path.join(model_dir, 'scaler.pkl'), 'rb') as f:
        scaler = pickle.load(f)
    data_scaled = pd.DataFrame(scaler.transform(data), index=data.index)
    return data_scaled


def ApplyTransformPCA(model_dir, data):
    with open(os.path.join(model_dir, 'pca.pkl'), 'rb') as f:
        pca = pickle.load(f)
    data_pca = pd.DataFrame(pca.transform(data), index=data.index)
    return data_pca


#def Predict(inhibitor, normed_specimen, pkl_1, pkl_2):
#    """Uses the pickled model to predict the AUC for the specimen."""
#    z_scores = (normed_specimen[pkl_1.gene] - pkl_1.gene_mean) / pkl_1.gene_std
#    return z_scores.dot(pkl_1[inhibitor]) + pkl_2.loc[inhibitor].intercept


def RunPredictions(model_dir, input_dir, output_dir):
    print('Loading models...')
    with open(os.path.join(model_dir, 'regressors.pkl'), 'rb') as f:
        regressors = pickle.load(f)

    print('Loading and pre-processing...')
    rnaseq = pd.read_csv(os.path.join(input_dir, 'rnaseq.csv'))

    specimens = TransposeRnaSeqTable(rnaseq)
    normed_specimens = NormSpecimens(specimens)
    selected = ApplySelectedGenes(model_dir, normed_specimens)
    scaled = ApplyStandardScaler(model_dir, selected)

    X = ApplyTransformPCA(model_dir, scaled)

    print('Predicting per-specimen AUC...')
    predictions = pd.DataFrame()

    for inhibitor in regressors.keys():
        regr = regressors[inhibitor]
        aucs = pd.DataFrame(regr.predict(X), index=X.index).reset_index()
        aucs.columns = ['lab_id', 'auc']
        aucs.insert(0, 'inhibitor', inhibitor)
        predictions = predictions.append(aucs)

    predictions.to_csv(os.path.join(output_dir, 'predictions.csv'), index=False)

    #    print('Getting the cartesian product of inhibitors and specimens...')
#    inhibitors = pkl_2.index
#    specimens = normed_specimens.index
#    aucs = pd.DataFrame(
#        product(inhibitors, specimens),
#        columns=['inhibitor', 'lab_id'])

#    print('Predicting per-specimen AUC...')
#    aucs['auc'] = aucs.apply(lambda r: (
#        Predict(r['inhibitor'], normed_specimens.loc[r['lab_id']], pkl_1, pkl_2)),
#                             axis=1)
