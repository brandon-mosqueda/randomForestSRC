# randomForestSRC

This repository is completely based on R's library [randomForestSRC](https://github.com/kogalur/randomForestSRC) [[1]](#1), it includes all the original source code with the custom splitting rules proposed and develop by [[2]](#2). Such new splitting rules are written in C and needs to be compiled with the whole project, for this reason we created this repository for automatically install the randomForestSRC library with these new splitting rules.

## Splitting rules

As described in [[2]](#2):

> We propose a general hurdle methodology to model a response from a homogeneous or a non-homogeneous Poisson process with excess zeros, based on two forests. The first forest in the two parts model is used to estimate the probability of having a zero. The second forest is used to estimate the Poisson parameter(s), using only the observations with at least one event.

> To build the trees in the second forest, we propose specialized splitting criteria derived from the zero truncated homogeneous and non-homogeneous Poisson likelihood. The particular case of a homogeneous process is investigated in details to stress out the advantages of the proposed method over the existing ones. Simulation studies show that the proposed methods perform well in hurdle (zero-altered) and zero-inflated settings.

## Our contribution

Provide a wrapper for implementing the Zero Altered Poisson (ZAP) random forest models proposed by [[2]](#2) for general purpose with count data with excess of zeros, a wrapper for tuning ZAP random forest models and a function of ZAP random forest models for genomic prediction.

## Authors

* Osval Antonio Montesinos López (Author)
* Abelardo Montesinos López (Author)
* Brandon A. Mosqueda Gonzalez (Author, Maintainer)
* José Cricelio Montesinos López (Author)
* José Crossa (Author)
* Carlos Alberto Flores Cortés (Author)

## References

<a id="1">[1]</a> Ishwaran H. and Kogalur U.B. (2020). Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC), R package version 2.9.3.

<a id="2">[2]</a> Mathlouthi, W., Larocque, D., & Fredette, M. (2020). Random forests for homogeneous and non-homogeneous Poisson processes with excess zeros. Statistical Methods in Medical Research, 29(8), 2217–2237. https://doi.org/10.1177/0962280219888741