<!-- README.md is generated from README.Rmd. Please edit that file -->


SPIEC-EASI
=========

Sparse InversE Covariance estimation for Ecological Association and Statistical Inference

This is the SPIEC-EASI MATLAB package which will be useful to anybody who wants to infer graphical models for absolute abundance and compositional data. Its primary application is intended for microbiome relative abundance data (generated from 16S amplicon sequence data). It also includes a generator for [overdispersed, zero inflated] multivariate, correlated count data. Please see the paper published in [PLoS Comp Bio](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). 

For R code please see https://github.com/zdk123/SpiecEasi

Since the original publications several additions have been made to the package. In particular, a lower bound on the lambda-path is estimated from subsamples for speeding up the stability selection part. This procedure is detailed in
[Generalized Stability Approach for Regularized Graphical Models](https://arxiv.org/abs/1605.07072).

Neighborhood selection can be performed both with the Lasso and (approximate) TREX. The latter approach is based on the methods detailed in 
[Topology Adaptive Graph Estimation in High Dimensions](https://arxiv.org/pdf/1410.7279.pdf)

## Installation ##

The package is self-contained. Code for fast L1-constraint linear regression is added from Mark Schmidt's homepage 
https://www.cs.ubc.ca/~schmidtm/Software/L1General.html under /solvers/schmidt-projected-subgradient. 

## Basic Usage ##

In the /examples/ folder you find several examples about the different modes of usage. 

## Analysis of a subset of American Gut data ##


