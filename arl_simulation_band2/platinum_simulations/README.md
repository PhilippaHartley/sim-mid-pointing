These are the platinum simulations.

From https://confluence.skatelescope.org/display/SE/SIM+Team%3A++Agree+pointing+simulation+parameters+%28SIM-43%29+next+steps+-+280619

Summary:
 - 1.4GHz, full track 45 degree declination, +/- 6 hours
 - All cases use GRASP beam at current FOV.
 - Wind heading tracked throughout simulation, using nearest (in az, el) PSD
 - Reference pointing at 30min

The variant cases are:

 - Case 1: Scaling with wind speeds. Wind heading tracked throughout simulation, using nearest (in az, el) PSD. 
 Frequency is L band. 0.3Jy/beam cutoff/ 10s integration
 - Case 2. Fixed offsets in elevation and cross-elevation. Pointing mis-set due to 
pointing error. Equivalent to case 1 but with same RMS but errors constant.
 - Case 3: Snapshot version of case 1
 - Case 4: There is no case 4
 - Case 5: Variation with integration time (30/3)
 - Case 6: Variation with elevation (i.e. change declination: +15/0/-85)
 - Case 7: Case 1 with MID_GAUSS