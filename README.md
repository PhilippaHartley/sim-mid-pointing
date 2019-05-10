# sim-mid-pointing

This repository contains three different sets of software for running SKA pointing error simulations.

## ARL

...to be added 

## Meqtrees

Requires the installation of meqtrees: http://meqtrees.net/ and CASA: https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing.

Requires Python (2 or 3)

#### Usage
```
cd meqtrees_simulation
python run_meqtrees.py
```

## CASA

Install CASA from https://casa.nrao.edu/casadocs-devel/stable/introduction/obtaining-and-installing. 

Requires Python (2 or 3)

#### Usage
```
cd casa_simulation
casa -c pointing_errorsSKA.py
```



