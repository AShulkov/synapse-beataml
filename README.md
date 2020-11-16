Here is the results of a brief research for CTD-squared BeatAML DREAM Challenge ([syn20940518](syn20940518)) 

## Summary

It was analyzed if the baseline solution could be improved by feature preprocessing using dimensionality reduction algorithms. 


## Introduction

Due to the limited number of data points in training data and an extremely large quantity of independent variables (213 specimen vs 63677 genes only in gene expression table), a problem of decreasing the number of acting features and overfitting prevention seems to be one of the most crucial.  In this light methods of dimensionality reduction could be promising among other feature selection algorithms. They are widely used to effectively squeeze feature space while saving most of the information сontained in the initial data. Along with сommon linear methods (such as PCA - principal component analysis), more sophisticated non-linear methods were introduced, e.g. kernel PCA, MDS, and t-SNE. 


## Methods

First, the number of most variably expressed genes (the genes with the largest variance on the training dataset) was experimentally tuned using leave one out cross-validation and Ridge regression algorithm. It turned out that more accurate tuning of the number of most expressed genes (20000 vs initial 1000) saves more information for subsequent processing.
Second, a number of dimensionality reduction algorithms were tested: PCA and truncated SVD algorithms as well as their modifications - kernel PCA with the distinct types of kernels (linear, polynomial, RBF, sigmoid, cosine), and other non-linear algorithms (MDS and t-SNE). The resulting dataset was used as input for the Ridge regression algorithm. It was found that applying most of the methods gives an increase both to the default regression metric (negative mean squared error) and metric specific for this subchallenge (averaged Spearman correlation of predictions for each drug).  The best results were obtained while using the kernel PCA algorithm with RBF kernels.

All methods were implemented using Python's Scikit-learn library implementation.


## Conclusion

It looks like both linear and kernel PCA can be used for effective dimensionality reduction for gene expression data. The resulting dataset allows easier and faster analysis with the usage of different ML algorithms.  


## The Model

One Ridge Regression model is trained for each inhibitor to predict AUC. The only input is gene expression (rnaseq.csv).

Specifics:

* The 20000 most variable genes are used for training
* The log2(cpm) values are normalized per-specimen
* The z-score is computed for each gene
* Ridge Regression is trained using leave-one-out cross-validation to predict AUC


## On-Disk Representation

The trained model is stored in four "pickles":

- selected_genes.pkl: contains the list of 20000 the most expressed genes
- scaler.pkl: contains fitted scaler model (z-score computation for each gene)
- pca.pkl: contains fitted parameters of dimensionality reduction algorithm (KernelPCA)
- regressors.pkl: contains fitted Ridge Regression models for each inhibitor

## References

Bernhard Schoelkopf, Alexander Smola, Klaus-Robert Mueller, Kernel principal component analysis, In Advances in kernel methods, MIT Press, Cambridge, MA, USA 327-352, 1999.

I. Borg, P. Groenen, Modern Multidimensional Scaling: Theory and Applications, Springer Series in Statistics, 1997.

L.J.P. van der Maaten, G.E. Hinton, Visualizing High-Dimensional Data Using t-SNE, Journal of Machine Learning Research 9:2579-2605, 2008.
