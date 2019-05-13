# sim-mid-pointing

This repository contains three different sets of software for running SKA pointing error simulations.

## ARL

...to be added 

## Meqtrees

Notes: The current output of the meqtrees script produces very similar results to the ARL output, except that the residuals look different. This could be an image weighting issue but needs checking.

Requires the installation of meqtrees: http://meqtrees.net/. Some libraries required for the installation are out of date; up-to-date libraries are included in the instructions for building from source on Ubuntu 18: https://confluence.skatelescope.org/display/SE/Installing+and+running+meqtrees

Requires the installation of CASA: https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing.



Requires Python (2 or 3)

#### Usage
```
cd meqtrees_simulation
python run_meqtrees.py
```

## CASA

Notes: CASA simulator tool contains a function for simulating pointing errors -- sm.setpointingerror -- but this function is not yet functional. The script below attempts a workround by modifiying the pointing table manually before predicting visibilities, but this does not work. The CASA EPJones table does not appear to be fully supported yet. 

Install CASA from https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing. 

Requires Python (2 or 3)

#### Usage
```
cd casa_simulation
casa -c pointing_errorsSKA.py
```



