These are the palladium simulations. Improved on platinum by fixing factor of two 
excess in pointing noise.

From https://confluence.skatelescope.org/display/SE/SIM+Team%3A++Agree+pointing+simulation+parameters+%28SIM-43%29+next+steps+-+280619

Summary:
 - 1.4GHz, full track 45 degree declination, +/- 6 hours
 - All cases use GRASP beam at current FOV.
 - Wind heading tracked throughout simulation, using nearest (in az, el) PSD
 - No Reference pointing

The main simulation is case5/10s. 

The variant cases are:

 - Case 2. Fixed offsets in elevation and cross-elevation. Equivalent to case5/10s 
 but with same RMS but errors constant.
 - Case 3: Snapshot version of case5/30s
 - Case 5: Variation with integration time (100/30/10/3)
 - Case 6: Variation with elevation (i.e. change declination: +15/0/-75)
 - Case 7: case5/30s with MID_GAUSS