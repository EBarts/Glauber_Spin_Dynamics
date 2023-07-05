# Glauber Spin Dynamics

This repository contains Monte Carlo methods for studying the non-equilibrium
behavior of spin systems. The codes provided here were used to study the magnetization
reversal process in Zn-doped kamiokites, as described in the paper
titled "Magnetization reversal through an antiferromagnetic state"
([arXiv link](https://arxiv.org/abs/2211.05028)). Theoretical results
presented in the paper can be reproduced using the data and plotting procedures in the "Data_and_Plotting_procedures" folder.

## Main Script

The main script allows us to calculate the averaged 
magnetization and electric polarization in a two-sublattice model of
Ising spins coupled to isotropic Heisenberg spins (see Eq.(3) in the
paper). To compile the code, use the following command:

```
g++ -o FswPar -O3 -fopenmp 3DquantumIsingFieldSweepPar.cpp
```

## Plotting

To visualize the observables, one can use the provided MATLAB script named
"plotWeightedAverages.m". Alternatively, after extracting the raw data
with the "saverawdata.m"
script in the "Data_and_Plotting_procedures" folder, one can use other software for visualization. Please
note that calculations were performed using the Peregrine
cluster at the University of Groningen using
"inputfile_Field_Sweep.sh" as
the input file.
Multiple field sweeps are collected in the "FSdata" folder, with the averaged observables being displayed in two PNG files.

## Contact

If you have any questions or require further clarification regarding
the procedure, feel free to reach out to the corresponding authors via
email. Additionally, the MATLAB version of the main script, as well as
other scripts utilized in my PhD thesis titled "Unconventional
Magnetic States and Defects" supervised by Prof. M. Mostovoy, will be uploaded to this repository soon.
