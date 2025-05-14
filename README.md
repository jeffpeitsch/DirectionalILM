# DirectionalILM

This repository contains code to simulate from and fit the directional ILMs described in "Directional Spatial Individual Level Models of Infectious Disease Transmission" by Peitsch, Pokharel, and Deardon; available at https://hdl.handle.net/1880/121536. 

SimStudyPopulation.txt provides the (X,Y) locations of individiduals used in the simulation study. 

SimulationVonMises.c and SimulationWrappedCauchy.c provide code to simulate epidemics using the directional ILM with their respective circular distributions. 

SimStudyVonMisesMCMC.c and SimStudyWrappedCauchyMCMC.c provide code to fit the directional ILM with their respective circular distributions. This code was used for the simulation study as well as the TSWV analyses, as the same version of directional ILM was used in both cases.

FMDVonMisesMCMC.c and FMDWrappedCauchyMCMC.c provide code to fit the directional FMD-ILM with their respective circular distributions. This code was used for the FMD analyses, as the directional FMD-ILM contains additional parameters specific to this dataset. 
