import numpy as np, os
from matplotlib import pyplot as plt
from scipy import optimize


psd = np.loadtxt("PSD_data/El45Az135.dat")

plotdir = "plots/"

if not os.path.isdir(plotdir):
    os.system('mkdir '+plotdir)

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

axes = ["az","el","pxel", "pel"]

for axis in axes:
    print ('*** Creating times series from PSD of %s axis ***'%axis)
    az = axesdict[axis]

    freq_interval = 0.001
    if (axis=="az") or (axis=="el"):   
        # determine index of maximum PSD value; add 50 for better fit
        az_max_index = np.argwhere(az==np.max(az))[0][0]+50
        max_freq = 0.4

        freq_max_index = np.argwhere(freq>max_freq)[0][0]
    else: 
        az_max_index = 400 # not max; just a break
        print ('frequency break: ', az_max_index)
        max_freq = 0.1
        freq_max_index = np.argwhere(freq>max_freq)[0][0]     
        print ('max frequency: ', max_freq)   

    # fit polynomial to psd up to max value
    p_az1 = np.polyfit(freq[:az_max_index], az[:az_max_index], 5)
    f_az1 = np.poly1d(p_az1)
    # fit polynomial to psd beyond max value
    p_az2 = np.polyfit(freq[az_max_index:freq_max_index], az[az_max_index:freq_max_index], 5)
    f_az2 = np.poly1d(p_az2)

    # resampled to construct regularly-spaced frequencies
    regular_freq1 = np.arange(freq[0], freq[az_max_index], freq_interval)
    regular_az1 = f_az1(regular_freq1)
    regular_freq2 = np.arange( freq[az_max_index],freq[freq_max_index] ,freq_interval)
    regular_az2 = f_az2(regular_freq2)
    regular_freq = np.append(regular_freq1, regular_freq2)
    regular_az = np.append(regular_az1,regular_az2)

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
    z_az = np.append(0+0*1j,z_az)
    regular_freq = np.append(0, regular_freq)

    # perform inverse fft
    ts = np.fft.ifft(z_az) 

    # set up and check scalings
    M=len(regular_freq)
    print ('Number of frequencies: ', M) # N=M+1

    N = len(ts)
    print ('number of times: ', N)

    Fmax = np.max(regular_freq)
    print ('Fmax: ',Fmax)

    Dt = 1/Fmax
    print ('Dt: ',Dt )

    tmax = (N-1)*Dt
    print('tmax: ', tmax)
    print ()

    ts.real*=N # the result is scaled by number of points in the signal, so multiply - real part - by this

    # The output of the iFFT will be a random time series on the finite 
    # (bounded, limited) time interval t = 0 to tmax = (N-1) X Dt, #
    # where Dt = 1 / (2 X Fmax)

    # scale to time interval
    times = np.arange(len(ts))*Dt

    # plot results
    plt.clf()
    plt.plot( times, np.real(ts), label = '%s axis'%axis)
    plt.xlabel('time, s')
    plt.ylabel('pointing error, as')
    plt.legend()
    plt.savefig(plotdir+'time_series_%s.png'%axis)
    
    # check fft returns original PSD (amplitudes plotted again)
    ts.real/=np.float(N)
    fft =  np.fft.fft(ts)
    plt.clf()
    plt.scatter(np.arange(len(fft)), np.absolute(fft), label = '%s axis'%axis) 
    plt.xlabel('mode')
    plt.ylabel('amplitude, as/SQRT(Hz)')
    plt.legend()
    plt.savefig(plotdir+'check_FFT_%s.png'%axis)
