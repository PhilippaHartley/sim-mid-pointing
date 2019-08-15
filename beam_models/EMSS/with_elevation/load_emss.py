import scipy.io
from matplotlib import pyplot as plt
import numpy as np, os
from astropy.io import fits
import matplotlib.colors as colors
import matplotlib.ticker as ticker     
import glob




def cart2pol(x,y):
    r = np.sqrt(x**2 +y**2)
    theta = np.arctan2(x,y)
    return(r,theta)

#filenames = glob.glob()

# grasp is 280 pix for 8 degrees diam

# have 120 degrees diam >> 120/8 * 280 = 4200

freq = 1360
beamsize = 4200 # even number

data_dir = '/Users/p.hartley/Dropbox (SKA)/pointing_errors/code/beam_conversion/beam_deformation/EMSS/SKADCBeamPatterns/2019_08_06_SKA_SPFB2/'
out_dir = 'out/'
plots_dir = 'plots/'

if not os.path.isdir(out_dir):
    os.system('mkdir '+out_dir)

if not os.path.isdir(plots_dir):
    os.system('mkdir '+plots_dir)    


elevations = [15, 45, 90]


for el in elevations:

    data_dir2 = '2019_08_06_SKA_%s_B2/'%el
    file_name = "B2_%s_%d.mat"%(el,freq)

    D = scipy.io.loadmat(data_dir+data_dir2+file_name)
    #print (D.items())
    th = D["th"].squeeze() # angle from bore sight with th=0 defining bore sight [deg]
    ph = D["ph"].squeeze() # angle in aperture plane with 0="up from vertex" and defines plane of symmetry [deg]
    JHH = D["Jqh"].squeeze() # Jones "voltage" pattern, arranged rows:th x columns:ph
    JVV = D["Jpv"].squeeze()
    JHV = D["Jqv"].squeeze()
    JVH = D["Jph"].squeeze()

    # keep theta in deg?
    print (th.shape,ph.shape)


    ph_rad = (ph/360.)*(2*np.pi)
    th_rad = (th/360.)*(2*np.pi)

    uflat = np.array([])
    vflat = np.array([])
    for counti, i in enumerate(th):
        for countj, j in enumerate(ph_rad):
            v = ((i)*np.sin(j))
            u = ((i)*np.cos(j))
            uflat = np.append(uflat, u)
            vflat = np.append(vflat,v)

    tharr, pharr = np.meshgrid((th), ph_rad) 

    x, y = np.meshgrid(np.linspace(np.min(vflat),np.max(vflat), np.int(beamsize)), np.linspace(np.min(uflat),np.max(uflat), beamsize))

    r, theta = cart2pol(x,y)

    azel = np.zeros((4,np.int(beamsize),np.int(beamsize)), dtype = complex)

    poldict =    {
          "JHH": JHH,
          "JVV": JVV,
          "JHV": JHV,
          "JVH": JVH
        }
        
    polarr = np.dstack((JHH, JVV,JHV,JVH))

    pols = ["JHH","JVV","JHV", "JVH"]

    for k , pol in enumerate(pols):
        for i in range(beamsize):
            for j in range(np.int(beamsize)):
       

                azel[k,j,i]= ((polarr[:,:,k][(np.abs(r[i,j]-th).argmin(),np.abs(theta[i,j]+np.pi-ph_rad).argmin())]))


    #flipped_azel = np.flip(azel, 1)             
    #azel = np.hstack((flipped_azel,azel))


    '''
    hdr = fits.Header()
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = '-32'
    hdr['NAXIS'] =
    hdr['NAXIS1'] =
    hdr['NAXIS2'] = 
    hdr['NAXIS3'] =
    hdr['NAXIS4'] =
    hdr['WCAXES'] =
    hdr['CRPIX1'] =
    hdr['CRPIX2'] =
    hdr['CRPIX3'] =
    hdr['CRPIX4'] =
    hdr['CDELT1'] =
    '''

    hdr = fits.Header.fromstring("""\
SIMPLE  =                    T / conforms to FITS standard                      
BITPIX  =                  -32 / array data type                                
NAXIS   =                    4 / number of array dimensions                     
NAXIS1  =                  %s                                                  
NAXIS2  =                  %s                                                  
NAXIS3  =                    4                                                  
NAXIS4  =                    1                                                  
WCSAXES =                    4 / Number of coordinate axes                      
CRPIX1  =                 %s / Pixel coordinate of reference point            
CRPIX2  =                 %s / Pixel coordinate of reference point            
CRPIX3  =                  1.0 / Pixel coordinate of reference point            
CRPIX4  =                  1.0 / Pixel coordinate of reference point            
CDELT1  =                %.3f / [deg] Coordinate increment at reference point  
CDELT2  =                 %.3f / [deg] Coordinate increment at reference point  
CDELT3  =                 -1.0 / Coordinate increment at reference point        
CDELT4  =         999999.99999 / [Hz] Coordinate increment at reference point   
CUNIT1  = 'deg'                / Units of coordinate increment and value        
CUNIT2  = 'deg'                / Units of coordinate increment and value        
CUNIT4  = 'Hz'                 / Units of coordinate increment and value        
CTYPE1  = 'AZELGEO long'       / Coordinate type code                           
CTYPE2  = 'AZELGEO lati'       / Coordinate type code                           
CTYPE3  = 'STOKES'             / Coordinate type code                           
CTYPE4  = 'FREQ'               / Frequency (linear)                             
CRVAL1  =                  0.0 / [deg] Coordinate value at reference point      
CRVAL2  =                  0.0 / [deg] Coordinate value at reference point      
CRVAL3  =                 -5.0 / Coordinate value at reference point            
CRVAL4  =         %s / [Hz] Coordinate value at reference point       
LONPOLE =                180.0 / [deg] Native longitude of celestial pole       
LATPOLE =                -35.0 / [deg] Native latitude of celestial pole        
END         
 """%(beamsize,beamsize,beamsize/2,(beamsize/2)+1, -(2*np.max(th)/beamsize), (2*np.max(th)/beamsize), (freq*1000)), sep='\n')
    
   # primary_hdu = fits.PrimaryHDU(header = hdr)

    hdu_new = fits.PrimaryHDU((np.real(azel)), header=hdr)
    hdu_new.writeto(out_dir+'real_%sMHz.fits'%file_name.strip('.mat'))

    hdu_new = fits.PrimaryHDU((np.imag(azel)))
    hdu_new.writeto(out_dir+'imag_%sMHz.fits'%file_name.strip('.mat'))

    hdu_new = fits.PrimaryHDU((np.absolute(azel)))
    hdu_new.writeto(out_dir+'amp_%sMHz.fits'%file_name.strip('.mat'))

    hdu_new = fits.PrimaryHDU((np.angle(azel)))
    hdu_new.writeto(out_dir+'phas_%sMHz.fits'%file_name.strip('.mat'))

    plt.clf()

    fig, axes = plt.subplots(nrows=1, ncols=4, sharex = True, sharey = True, figsize = (6,4))
    for count, ax in enumerate(axes.flat):
        polarrabs = np.absolute(polarr)
        polarrabs[:,:,count] = 10*np.log((polarrabs[:,:,count])/np.max(polarrabs[:,:,count])) # power
        im = ax.imshow((polarrabs[:,:,count]), vmin = -70, vmax = 0, origin = 'lower', extent = [np.min(ph),np.max(ph), np.min(th), np.max(th)], aspect = ((np.max(ph)/np.max(th))*len(th))/len(ph))
        #if (count==2) or (count==3):
        ax.set_xlabel(r'$\phi$ [deg]')
        if (count==0):# or (count==2):
            ax.set_ylabel(r'$\theta$ [deg]')
        ax.set_title(pols[count])    

    plt.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), use_gridspec=True)
    #plt.show()
    plt.savefig(plots_dir+'pol_amp_%sMHz.png'%(file_name))



    fig, axes = plt.subplots(nrows=1, ncols=4, sharex = True, sharey = True, figsize = (6,4))
    for count, ax in enumerate(axes.flat):
        amp_azel = np.angle(azel)
        # db = 10*np.log((amp_azel)/np.max(amp_azel)) # power dB not field
        im = ax.imshow(np.angle(polarr[:,:,count]), vmin = -np.pi, vmax = np.pi, origin = 'lower', extent = [np.min(ph),np.max(ph), np.min(th), np.max(th)], aspect = ((np.max(ph)/np.max(th))*len(th))/len(ph))
        #if (count==2) or (count==3):
        ax.set_xlabel(r'$\phi$ [deg]')
        if (count==0):# or (count==2):
            ax.set_ylabel(r'$\theta$ [deg]')
        ax.set_title(pols[count])    

    plt.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), use_gridspec=True)

    plt.savefig(plots_dir+'pol_phas_%sMHz.png'%(file_name))




    azel_types = [np.absolute(azel), np.imag(azel),np.real(azel) ,np.angle(azel)]
    azel_names = ['amp', 'imag', 'real', 'phase']

    for i in range(4):
        azel_type = azel_types[i]
        plt.clf()
        fig, axes = plt.subplots(nrows=2, ncols=2, sharex = True, sharey = True, figsize = (8,6))
        for count, ax in enumerate(axes.flat):
            if i==0:
                azel_type[count,:,:] = 10*np.log((azel_type[count,:,:])/np.max(azel_type[count,:,:])) # power not voltage
                vmin = -70
                vmax = 0
            else:
                vmax  = np.max(azel_type[0,:,:])    
                vmin  = np.min(azel_type[0,:,:])     
            norm = None
            if (i==1) or (i==2):
                norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,   vmin=-1.0, vmax=1.0)
            im = ax.imshow((azel_type[count,:,:]), vmin = vmin, vmax = vmax, origin = 'lower', extent = [-(np.max(th)), np.max(th), -(np.max(th)), np.max(th)], norm=norm)
            if (count==2) or (count==3):
                ax.set_xlabel(r'$\rho$cos($\phi$) [deg]')   
            if (count==0) or (count==2):
                ax.set_ylabel(r'$\rho$sin($\phi$) [deg]')
            ax.set_title(pols[count])  
        plt.tight_layout()
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), use_gridspec=True)
        if (i==1) or (i==2):
            pass
           # cbar.set_ticks(ticker.LogLocator(), update_ticks=True)
        cbar.ax.tick_params(size=0)    
        plt.savefig(plots_dir+'azel_%s_%sMHz.png'%(azel_names[i],file_name))






















        #sin theta is radius
        # phi is angle around 
