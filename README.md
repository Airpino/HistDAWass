# HistDAWass
An R package for histogram data analysis
Package: HistDAWass

Maintainer: Antonio Irpino <antonio.irpino@unina2.it>
Description: The package collects classes and methods for data analysis of
    histogram-valued data, i.e., data described by univariate histograms. The
    methods and the basic statistics for histogram-valued data are mainly based
    on the L2 Wasserstein metric between distributions, i.e., a Euclidean
    metric between quantile functions.
License: GPL (>=2)
List of methods
- Basic statistics for histogrm variables:
    - mean
    - standard deviation
    - covariance
    - correlation
- Data analysis and visualization techniques
    - Unsupervised classification
      - dynamic clustering, k-means
      - dynamic clustering using adaptive distances
      - Batch Kohonen SOM maps
      - Batch Kohonen SOM maps with adaptive distances
      - fuzzy c-means
      - fuzzy c-means with adaptive distances
      - Principal Component Analysis for one histogram variable
      - Principal Component Analysis for multiple histogram variables
- Regression models
    - Least squares estimation of a two components liner model (IRINO-VERDE regression)
- Time series analysis
    - Moving average smoothing of histogram time series (HTS)
    - Exponential smoothing of HTS
    - K-NN prediction of HTS
