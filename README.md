# cuda-puff

This project was done in the Cardiac Systems Pharmacology Lab of Dr. Eric Sobie, Icahn School of Medicine at Mount Sinai (New York, NY), and with the guidance of M.D./Ph.D. candidate Deanalisa Jones.

A calcium puff is a local release of calcium ions from a store such as the endoplasmic reticulum, or the sarcoplasmic reticulum in muscle cells. 
Calcium is released through clusters of inositol 1,4,5-trisphosphate receptors (IP3Rs). This intracellular calcium signalling is essential for processes such as contraction in muscle cells and fertilization of oocytes.

[Siekmann et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/22947927), created a Markov model of the kinetics of IP3Rs, taking their state-switching into account. In 2013, [Cao et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3852038/) created a new IP3R model incorporating important dynamic properties of the receptors, which the Siekmann model failed to account for. The widely accepted Cao model uses the Gillespie algorithm with a variable time step.

The Cao model was written in MATLAB. This project aims to drastically increase the speed of the simulation using CUDA, a parallel computing platform that performs general purpose computations using a graphics processing unit. This enables us to run many trials simulataneously and in a significantly shorter period of time.
The program runs 5 trials in a few minutes, compared to the original MATLAB version, which takes several hours to run one trial.

The CUDA C++ (.cu) file provided outputs data of interest in the form of CSV files which can then be read into MATLAB. This repository also includes a MATLAB script that generates several figures from the simulation data, including: the fluorescence time course, and histograms displaying the distribution of puff amplitudes, durations, and inter-puff intervals.

Members of the Sobie Lab will use this program to vary parameters / do parameter sensitivity analysis to study the underlying differences between calcium sparks and puffs. 
