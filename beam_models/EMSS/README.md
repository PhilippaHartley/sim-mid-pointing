EM analysis is being performed by SARAO, who have provided an initial set of far field patterns in Matlab files e.g. 
 `SKA_B2_965.mat`. A script, `import_beam_B1.py`, loads this MATLAB file which contains voltage patterns for four
  polarisations: JHH, JVV, JHV and JVH. Each voltage pattern is provided as a theta x phi array, where theta is the angle from the bore sight, with theta = 0 defining the bore sight [deg], and phi is the angle in aperture plane, with 0 = "up from vertex", and defines plane of symmetry [deg]. The script converts each voltage pattern into a rho\*sin phi x rho\*cos phi array, where rho = theta and therefore retains the value of radial sky angle from bore sight.

The voltage patterns are processed as follows:
 - The absolute peak in the voltage pattern is normalised to 1.0
 - The phase gradient of the beam pattern in the focal plane is fit and removed
 - The peak of the voltage pattern is shifted (by integral pixels) to be at the nominal centre.
 
 The beam patterns are interpolated in elevation in steps of 1 degree, using a spline. 