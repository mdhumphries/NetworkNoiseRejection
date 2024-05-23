Spectral estimation for detecting network structure
============================

MATLAB toolbox for finding low-dimensional structure in networks using spectral estimation

Accompanies the paper "Spectral estimation for detecting low-dimensional structure in networks using arbitrary null models" available at [PLoS One](https://doi.org/10.1371/journal.pone.0254057) and [arXiv](https://arxiv.org/abs/1901.04747)

A Python port of this toolbox by Thomas Delaney is here: https://github.com/thomasjdelaney/Network_Noise_Rejection_Python

The main tools here are:
(1) spectral estimation algorithm to find a network's low-dimensional subspace, if one exists
(2) node rejection in that subspace, to detect network nodes that do not contribute to the low-dimensional network description
(3) consensus clustering to find network communities/modules in that subspace

How to apply this code to a NxN correlation matrix or other matrices of pairwise comparisons to study e.g. time-series correlations:
- we interpret the correlation values as weights of links in the network between the N nodes: so any comparison measure must be a measure of similarity between the N items
- the algorithm assumes network weights are non-negative: typically we rectify negative correlations to zero.
- the algorithm assumes no self-links, so diagonal weights/correlations should be zero.

## Code
The core functions are in the folder Network_Spectra_Functions/
The top-level script *Toolbox_Examples_Script* illustrates calling these functions on an example network, and visualising some of the results

Much of the rest of the code is:
* Community detection algorithms (Network_Analysis_Functions/ and ZhangNewman2015/)
* Code for generating synthetic networks with modular structure and noise (SyntheticModel/)
* Scripts for generating all analyses in the above paper, and some preliminary visualisations (Scripts_For_Paper/)
* The set of data networks analysed in the above paper, including code for their generation where necessary (Networks/)
* The results of our analyses used in the paper (Results/)

### Features in the Toolbox that are not in the paper
The toolbox was built as a general purpose code base for exploring the potential of spectral rejection. Thus the main functions have a number of options that have yet to be explored. These include:

#### The lower bound on the predicted eigenvalues from the null model
In the paper, we examine network structure found when eigenvalues in the data network exceed the upper bound predicted by the null model. These eigenvalues describe low-dimensional clusters. However, we could equally ask if eigenvalues in the data network fall below (are more negative than) the lower bound predicted by the null model.

This option is available in _LowDSpace_

We demonstrate its use in the script *Negative_Eigenvalues* (in Scripts_For_Paper/)

#### Estimating the maximum eigenvalue
In the paper, we use the expectation over all sampled null models as our estimate of the maximum and minimum predicted eigenvalues. The function _LowDSpace_ (and so also _NodeRejection_) contains options for instead specifying an upper confidence interval on that maximum, computed over all sampled models. It also contains an option for [non-parametric prediction intervals](https://en.wikipedia.org/wiki/Prediction_interval). 

In the Supplementary Material, as a suggestion for future work we describe briefly three approaches to statistical testing of the data network eigenvalues, given the sample of null model maxium eigenvalues. Demonstrations of these are in the script *demo_of_statistical_tests* (in Scripts_For_Paper/)

##### Node projection
In the paper, we project nodes into a low-dimensional space using their L2-norm weighted by the eigenvalue of each dimension. We also test for node rejection using the expectation of the projections over all sampled null models. Function _NodeRejection_ has options for different weightings of the projections, and also for setting confidence and prediction intervals, as above.

