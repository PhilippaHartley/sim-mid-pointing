

Power spectra are derived from precision, standard and degraded wind conditions (< 5, <7 and <10 m/s, respectively). Directory ```PSD_data/``` contains sets of power spectra for each condition, for a selection of dish elevation angles and wind angles of attack. The filename pattern is El##Az##.dat, where El is the elevation abgle and Az is the wind angle of atcack. Within each set, axes EL and AZ refer to trajectory analysis of errors due to the servo movement in elevation and azimuth, respectively. PXEL and PEL refer to dynamic wind analysis errors in azimuth and elevation, respectively, where PXEL is the cross elevation in the azimuth direction, corrected for cos(EL).

PSD are provided for three elevation angles: EL15, EL45 and EL90. Calculations are performed for five directions of wind: AZ0 AZ45 AZ90 AZ135 AZ180, starting with AZ0 which means that the wind is blowing directly into the primary.

A script, time_series.py, generates a time series of pointing errors from a given power spectrum. 

From MTM: The simulation done in a way that the dynamic wind analysis were performed completely in the frequency domain. In contrast, trajectory simulations were done without wind. The reason for this was an intention to see the influence of trajectory dynamics (e.g. settling time) on the controller-structure couple separately from the influence of the gust wind.
