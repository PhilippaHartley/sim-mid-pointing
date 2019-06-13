import numpy as np, os, sys
from matplotlib import pyplot as plt
from scipy import optimize



psd = np.loadtxt("PSD_data/El45Az135.dat")

plotdir = "plots/"
outdir = "out/"

if not os.path.isdir(plotdir):
    os.system('mkdir '+plotdir)

if not os.path.isdir(outdir):
    os.system('mkdir '+outdir)

# inspect psds
plt.subplot(221)
plt.semilogx(psd[:,0], psd[:,1],label = 'S_az(f)')
plt.ylabel('as^2/Hz')
plt.legend()
plt.subplot(222)
plt.semilogx(psd[:,0], psd[:,2],label = 'S_el(f)')
plt.legend()
plt.subplot(223)
plt.semilogx(psd[:,0], psd[:,3], label = 'S_pxel(f)')
plt.ylabel('as^2/Hz')
plt.xlabel('f, Hz')
plt.legend()
plt.subplot(224)
plt.semilogx(psd[:,0], psd[:,4],label = 'S_pel(f)')
plt.legend()
plt.xlabel('f, Hz')
plt.savefig(plotdir+'power_spectra.png')

# define some arrays
freq = psd[:,0]
axesdict =    {
  "az": psd[:,1],
  "el": psd[:,2],
  "pxel": psd[:,3],
  "pel": psd[:,4]
}

doplot = 1 # turn off for speed

axes = ["az","el","pxel", "pel"]


time_stamps = 65 # currently need to increase (reduce) freq_interval to reduce (increase) number of time series points

for axis in axes:
    print ('*** Creating times series from PSD of %s axis ***'%axis)
    az = axesdict[axis]

    freq_interval = 0.001

    # the orginal PSD need to be resampled in order to obtain fixed frequency intervals fot the IFFT
    # each PSD is broken into two parts in order to produce better polynomial fits to the data
    # the az and el PSDs have a peak value mid-spectrum; the break is therefore located near this peak (trial and error to find best break)
    # the pxel and pel PSDs do not have a peak mid-spectrum, so a break is chosen at an arbitrary freequency after trial and error
    if (axis=="az") or (axis=="el"):   
        # determine index of maximum PSD value; add 50 for better fit
        az_max_index = np.argwhere(az==np.max(az))[0][0]+50
        # the maximum frequency is chosen to represent the top of the range of significant power; this may need to change for different PSDs
        # at low power values there is a danger of the fit producing negative values at the re-sampling stage
        max_freq = 0.4
        freq_max_index = np.argwhere(freq>max_freq)[0][0]
    else: 
        break_freq = 0.01 # not max; just a break
        az_max_index = np.argwhere(freq>break_freq)[0][0]
        # the maximum frequency is chosen to represent the top of the range of significant power; this may need to change for different PSDs
        # at low power values there is a danger of the fit producing negative values at the re-sampling stage
        max_freq = 0.1
        freq_max_index = np.argwhere(freq>max_freq)[0][0]     

    print ('Frequency break: ', freq[az_max_index])
    print ('Max frequency: ', max_freq)  

    if az_max_index>=freq_max_index:
        print ('Frequency break is higher than highest frequency; select a lower break')    
        sys.exit()

    # fit polynomial to psd up to max value
    p_az1 = np.polyfit(freq[:az_max_index], az[:az_max_index], 5)
    f_az1 = np.poly1d(p_az1)
    # fit polynomial to psd beyond max value
    p_az2 = np.polyfit(freq[az_max_index:freq_max_index], az[az_max_index:freq_max_index], 5)
    f_az2 = np.poly1d(p_az2)


    #freq_interval2 = (freq[freq_max_index]-freq[0])/time_stamps

    # resampled to construct regularly-spaced frequencies
    regular_freq1 = np.arange(freq[0], freq[az_max_index], freq_interval)
    regular_az1 = f_az1(regular_freq1)
    regular_freq2 = np.arange( freq[az_max_index],freq[freq_max_index] ,freq_interval)
    regular_az2 = f_az2(regular_freq2)
    regular_freq = np.append(regular_freq1, regular_freq2)
    regular_az = np.append(regular_az1,regular_az2)

    if (regular_az<0).any():
        print ('Resampling returns negative power values; change fit range')
        sys.exit()

    if doplot:
        # check sampling matches psd
        plt.clf()
        plt.semilogx(freq, (az), label = 'same')
        plt.semilogx(regular_freq, regular_az, label = 'regular')
        plt.ylabel('PSD, as^2/Hz')
        plt.xlabel('f, Hz')
        plt.legend()    
        plt.savefig(plotdir+'sampled_psd_%s.png'%axis)

    # get amplitudes from psd values
    amp_az = np.sqrt(regular_az) 

    # window function to smooth edges
    #amp_az = np.hamming(len(amp_az))*amp_az # not sure if necessary?

    # generate some random phases
    phi_az = np.random.rand(len(regular_az))*2*np.pi
  
    # create complex array
    z_az = amp_az * np.exp(1j*phi_az) # polar

    # construct Hermitian PSD in order to produce real signal
 
    # make symmetrical frequencies
    mirror_z_az = np.copy(z_az)

    # make complex conjugates
    mirror_z_az.imag-=2*z_az.imag

    # make negative frequencies
    mirror_regular_freq = -regular_freq 

    # join
    z_az = np.append(z_az,mirror_z_az[::-1])
    regular_freq = np.append(regular_freq,mirror_regular_freq[::-1])

    # add a 0 Fourier term
    z_az = np.append(0.* np.exp(1j*0.),z_az) # 0 DC component since mean 0 of time series(?)
    regular_freq = np.append(0, regular_freq)
    # set up and check scalings
    M=len(z_az)
    print ('Number of frequencies: ', M) # N=M+1


    # perform inverse fft
    ts = np.fft.ifft(z_az) 


    N = len(ts)
    print ('number of times: ', N)


    Fmax = np.max(regular_freq)
    print ('Fmax: ',Fmax)

    Dt = 1/Fmax
    print ('Dt: ',Dt )
#
    tmax = (N-1)*Dt
    print('tmax: ', tmax)
    print ()

    ts.real*=np.sqrt(N) # the result is scaled by number of points in the signal, so multiply - real part - by this
  # see https://docs.scipy.org/doc/numpy/reference/routines.fft.html for scaling desciption
  # ! this is wrong - appears to be scaled by 1/sqrt(N) - maybe since hermitian?


    # The output of the iFFT will be a random time series on the finite 
    # (bounded, limited) time interval t = 0 to tmax = (N-1) X Dt, #
    # where Dt = 1 / (Fmax)

    # scale to time interval
    times = np.arange(len(ts))*Dt

    if doplot:

        # plot results
        plt.clf()
        plt.plot( times, np.real(ts), label = '%s axis'%axis)
        plt.xlabel('time, s')
        plt.ylabel('pointing error, as')
        plt.legend()
        plt.title(  'r.m.s. = %s'%(np.std(np.real(ts))))
        plt.savefig(plotdir+'time_series_%s.png'%axis)
    
    np.save(outdir+'time_series_%s'%axis, ts.real)

    if doplot:
        # check fft returns original PSD (amplitudes plotted again)
        ts.real/=np.sqrt(N)
        fft =  np.fft.fft(ts)

        plt.clf()
        plt.scatter(np.arange(len(fft)), (np.absolute(fft))**2, label = '%s axis'%axis) 
        plt.xlabel('mode')
        plt.ylabel('amplitude, as/SQRT(Hz)')
        plt.legend()
        plt.savefig(plotdir+'check_FFT_%s.png'%axis)
