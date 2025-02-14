# Motor Cortex Simultaneous Coding
 MATLAB code that generates the figures from the paper
 
> Motor Cortex Latent Dynamics Encode Spatial and Temporal Arm Movement Parameters Independently
> 
> Andrea Colins Rodriguez, Matt G. Perich, Lee E. Miller, Mark D. Humphries
> 
> Journal of Neuroscience 28 August 2024, 44 (35) e1777232024; [DOI: 10.1523/JNEUROSCI.1777-23.2024](https://www.jneurosci.org/content/44/35/e1777232024)

## Summary

The main objective of this code is to perform PCA on the neural activity of PMd and M1 recorded while non-human primates performed a sequential movement task (see [Lawlor et al,. 2018](https://doi.org/10.1007/s10827-018-0696-6) for a full description of the task). 
Movements and neural activity are binned into 32 conditions (4 duration bins and 8 direction bins). 

## Dataset

All recordings from the Random-Target task can be found on [DANDI](https://dandiarchive.org/dandiset/000688/draft). The recordings used for the analyses are the following:
|                **Monkey M**                |
|:------------------------------------------:|
| sub-M_ses-RT-20140116_behavior+ecephys.nwb |
| sub-M_ses-RT-20140214_behavior+ecephys.nwb |
| sub-M_ses-RT-20140221_behavior+ecephys.nwb |
|                **Monkey C**                |
| sub-C_ses-RT-20150318_behavior+ecephys.nwb |
| sub-C_ses-RT-20150320_behavior+ecephys.nwb |
| sub-C_ses-RT-20150317_behavior+ecephys.nwb |
| sub-C_ses-RT-20150316_behavior+ecephys.nwb |
|                **Monkey T**                |
| sub-T_ses-RT-20130904_behavior+ecephys.nwb |
| sub-T_ses-RT-20130906_behavior+ecephys.nwb |
| sub-T_ses-RT-20130910_behavior+ecephys.nwb |

## Running the analysis

The script Main_analyses.m sets all the initial parameters and executes the code to perform PCA on all the neural recordings

1) Download the all the recording described in the table above and the code of this repository
2) Install the Matlab interface to read data from Neurodata Without Borders ([MATNWB](https://github.com/NeurodataWithoutBorders/matnwb))
3) Open the file reading_data_from_DANDI.m and set the path where MATNWB is stored.
4) Download the [dPCA toolbox from its repository](https://github.com/machenslab/dPCA)
5) Set the Matlab Path to the folder Data
6) Open Main_analyses.m and define the paths where the dPCA toolbox and this code are
7) Run adding to the current path.

## Expected output

This script generates:

- 6 figures corresponding to the main figures of the text.

- supplementary figures that were used to compute complementary statistics.
