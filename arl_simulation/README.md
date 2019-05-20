
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