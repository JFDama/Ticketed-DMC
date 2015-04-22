# Ticketed-DMC
A generic educational implementation of the Ticketed Diffusion Monte Carlo algorithm in C++

[![DOI](https://zenodo.org/badge/9297/JFDama/Ticketed-DMC.svg)](http://dx.doi.org/10.5281/zenodo.17001)

---
#### Algorithm

[Ticketed Diffusion Monte Carlo (TDMC)](http://arxiv.org/abs/1207.2866) is a recent Monte Carlo sampling alogrithm conceived by [Martin Hairer](http://www.hairer.org/) and [Jonathan Weare](http://www.stat.uchicago.edu/~weare/). It estimates the same types of quantities that Diffusion Monte Carlo does, but it does so with lower variance in every case, especially for small timesteps and with birth/death weight functions that change quickly. The difference is so striking in those cases that TDMC has new application areas that would barely be imaginable with vanilla DMC.

Most importantly for Hairer and Weare's paper, it looks remarkably suitable for calculating finite-time expectations of rare events and for statistical data assimilation.

---
#### This Code

This repo provides a C++11 implementation of a generic TDMC algorithm with the type of state to sample, the sampling dynamics, and the birth/death weight function as template parameters. It's designed to be easy to play around with above all; the performance is moderate and the implementation isn't as cache-friendly as it could be for the sake of code clarity. It's one header file, half comments, that requires only standard libraries.

Along with the algorithm code, the repo also contains a copy of the Hairer and Weare ArXiv paper introducing the algorithm, a correction to the paper, and a couple of simple examples including the LJ7 cluster example from the Hairer and Weare manuscript. Please try the examples, compare your results to the (corrected) reference values, and experiment!
