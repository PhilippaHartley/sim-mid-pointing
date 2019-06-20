

The power spectra are derived from precision wind conditions (< 5 m/s)Directory ```PSD_data/precision``` contains sets of PSDs for a selection of dish elevation angles and wind angles of attack. The filename pattern is El##Az##.dat, where El is the elevation abgle and Az is the wind angle of atcack. Within each set, axes EL and AZ refer to trajectory analysis of errors due to the servo movement in elevation and azimuth, respectively. PXEL and PEL refer to dynamic wind analysis errors in azimuth and elevation, respectively, where PXEL is the cross elevation in the azimuth direction, corrected for cos(EL).


A script, time_series.py, generates a time series of pointing errors from a given power spectrum. 


From MTM: The simulation done in a way that the dynamic wind analysis were performed completely in the frequency domain. In contrast, trajectory simulations were done without wind. The reason for this was an intention to see the influence of trajectory dynamics (e.g. settling time) on the controller-structure couple separately from the influence of the gust wind.
