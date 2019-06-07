# sim-mid-pointing

This repository contains three different sets of software for running SKA pointing error simulations.

## ARL

Requires the installation of the Algorithm Reference Library: https://github.com/SKA-ScienceDataProcessor/algorithm-reference-library. A step-through of the installation process is available at https://confluence.skatelescope.org/display/SE/Installing+and+Running+ARL.


## Meqtrees

Requires the installation of meqtrees: http://meqtrees.net/. Some libraries required for the installation are out of date; up-to-date libraries are included in the instructions for building from source on Ubuntu 18: https://confluence.skatelescope.org/display/SE/Installing+and+running+meqtrees

Requires the installation of CASA: https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing.


## CASA

Notes: CASA simulator tool contains a function for simulating pointing errors -- sm.setpointingerror -- but this function is not yet functional. The script below attempts a workround by modifiying the pointing table manually before predicting visibilities, but this does not work. The CASA EPJones table does not appear to be fully supported yet. 

Install CASA from https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing. 








