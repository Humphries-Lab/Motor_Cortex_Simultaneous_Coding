# Motor_Cortex_Simultaneous_Coding
 MATLAB code that generates the figures from the paper "Independent latent encoding of arm movement direction and extent in motor cortex"
 

## Summary

The main objective of this code is to perform PCA on the neural activity of PMd and M1 recorded while non-human primates performed a sequential movement task (see [Lawlor et al,. 2018](https://doi.org/10.1007/s10827-018-0696-6). 'Population coding of conditional probability distributions in dorsal premotor cortex' for a full description of the task). Movements and neural activity are binned into 32 conditions (4 duration bins and 8 direction bins). 
A full description of the datset format can be found on [here](http://crcns.org/data-sets/motor-cortex/pmd-1)


## Running the analysis

The script Main_analyses.m sets all the initial parameters and executes the code to perform PCA on all the neural recordings

1) Download the folders Data and Code
2) Set the Matlab Path to the folder Data
3) Open Main_analyses.m and run adding to the path.

## Expected output

This script generates 13 figures:

7 figures from the main text (Figure 5 of the paper is divided into 2 Matlab figures)
5 supplementary figures.
