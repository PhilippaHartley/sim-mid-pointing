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