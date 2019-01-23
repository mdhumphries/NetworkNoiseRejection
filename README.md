Network Noise Rejection
============================

MATLAB toolbox for finding and rejecting noise in networks

Accompanies the paper "Spectral rejection for testing hypotheses of structure in networks": https://arxiv.org/abs/1901.04747

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

##### Node projection
In the paper, we project nodes into a low-dimensional space using their L2-norm weighted by the eigenvalue of each dimension. We also test for node rejection using the expectation of the projections over all sampled null models. Function _NodeRejection_ has options for different weightings of the projections, and also for setting confidence and prediction intervals, as above.

