

The power spectra are derived from precision wind conditions (< 5 m/s), measured for an elevation and azimuth of 45 deg and 135 deg, respectively (the values represent the worst cases for axes separately i.e. do not correspond to the same position of the dish).


A script to generate a time series of pointing errors from a given power spectrum. Axes EL and AZ refer trajectory analysis of errors due to the servo movement in elevation and azimuth, respectively. PXEL and PEL refer to dynamic wind analysis errors in azimuth and elevation, respectively, where PXEL is the cross elevation in the azimuth direction, corrected for cos(EL).

The power spectra are derived from precision wind conditions (< 5 m/s), measured for an elevation and azimuth of 45 deg and 135 deg, respectively (the values represent the worst cases for axes separately i.e. do not correspond to the same position of the dish).

From MTM: The simulation done in a way that the dynamic wind analysis were performed completely in the frequency domain. In contrast, trajectory simulations were done without wind. The reason for this was an intention to see the influence of trajectory dynamics (e.g. settling time) on the controller-structure couple separately from the influence of the gust wind.
