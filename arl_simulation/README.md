**ARL Simulations of the effect of pointing observations on mid observations**

An observation with MID is simulated:

  - +/- 6 hours observation of a source at declination-45 deg
  - Frequency 1.4GHz 
  - pointing errors per antenna per integration added
  - Beam model is an unblocked, uniformly illuminated disk
  
To run an instance, make a directory e.g. simulation97 cp do_sim.sh, change 
the seed if appropriate and execute it using:

    sh dosim.sh
    
This loops through pointing errors of 1, 2, 4..., 128, 256 arcsec. A residual image
is produced for each, and written as png and fits files.

The statistics of each image are gathered and written into a csv file and plotted 
as well.

Each directory contains a driver script using command line parameters to control the sim. The same driver program, 
pointing_simulation, works for all case. It is copied into the simulation directory at the initialisation.

sims 0, 1, 2, 3: -6 to + 6 hours at ha=15 deg, dec=45 deg, different seeds
sim 4: same as sim 0 but snapshot
sim 5: same as sim 0 but static errors
sim 6: same as sim 0 but source on opposite side of pointing centre
sim 7: s3sky down to 0.3Jy with MID_GAUSS beam
sim 8: s3sky down to 0.1Jy with MID_GAUSS beam
sim 9: s3sky down to 0.03Jy with MID_GAUSS beam
sim 10: s3sky down to 0.03Jy with MID_GAUSS beam and static errors
sim 11: repeat of sim 0 with calculations in AZELGEO
sim 12: repeat of sim 5 with calculations in AZELGEO
sim 13: repeat of sim 0 with GRASP beam and calculations in AZELGEO
sim 14: repeat of 0 with az offset multiplied by cos(dec) - to mimic MeqTrees
sim 15: repeat of 0 with az offsets set to zero
sim 16: repeat of 0 with el offsets set to zero
sim 17: repeat of 0 in RADEC frame (use_radec=True)