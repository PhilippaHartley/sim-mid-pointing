import scipy.io
from matplotlib import pyplot as plt
import numpy as np, os
from astropy.io import fits
import matplotlib.colors as colors
import matplotlib.ticker as ticker     
def cart2pol(x,y):
    r = np.sqrt(x**2 +y**2)
    theta = np.arctan2(x,y)
    return(r,theta)


freq = 965
beamsize = 512 # even number


D = scipy.io.loadmat("SKA_B2_%d.mat"%freq)
th = D["th"].squeeze() # angle from bore sight with th=0 defining bore sight [deg]
ph = D["ph"].squeeze() # angle in aperture plane with 0="up from vertex" and defines plane of symmetry [deg]
JHH = D["Jqh"].squeeze() # Jones "voltage" pattern, arranged rows:th x columns:ph
JVV = D["Jpv"].squeeze()
JHV = D["Jqv"].squeeze()
JVH = D["Jph"].squeeze()

# keep theta in deg?
print (ph.shape,th.shape)
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

x, y = np.meshgrid(np.linspace(np.min(vflat),np.max(vflat), np.int(beamsize/2)), np.linspace(np.min(uflat),np.max(uflat), beamsize))
r, theta = cart2pol(x,y)



azel = np.zeros((4,np.int(beamsize/2),np.int(beamsize)), dtype = complex)

poldict =    {
      "JHH": JHH,
      "JVV": JVV,
      "JHV": JHV,
      "JVH": JVH
    }
    
    
polarr = np.dstack((JHH, JVV,JHV,JVH))
    

pols = ["JHH","JVV","JHV", "JVH"]

for k , pol in enumerate(pols):
    print (poldict[pol].shape)
    print (polarr[:,:,k].shape)
    for i in range(beamsize):
        for j in range(np.int(beamsize/2)):

            azel[k,j,i]= ((polarr[:,:,k][(np.abs(r[i,j]-th).argmin(),np.abs(theta[i,j]-ph_rad).argmin())]))


flipped_azel = np.flip(azel, 1)             
azel = np.hstack((flipped_azel,azel))

hdu_new = fits.PrimaryHDU((np.real(azel)))
hdu_new.writeto('real_%s.fits'%freq)

hdu_new = fits.PrimaryHDU((np.imag(azel)))
hdu_new.writeto('imag_%s.fits'%freq)

hdu_new = fits.PrimaryHDU((np.absolute(azel)))
hdu_new.writeto('amp_%s.fits'%freq)

hdu_new = fits.PrimaryHDU((np.angle(azel)))
hdu_new.writeto('phas_%s.fits'%freq)

plt.clf()

fig, axes = plt.subplots(nrows=1, ncols=4, sharex = True, sharey = True, figsize = (6,4))
for count, ax in enumerate(axes.flat):
    polarrabs = np.absolute(polarr)
    polarrabs[:,:,count] = 10*np.log((polarrabs[:,:,count])/np.max(polarrabs[:,:,count])) # pow
    im = ax.imshow((polarrabs[:,:,count]), vmin = -70, vmax = 0, origin = 'lower', extent = [0,180, 0, 6], aspect = (30*501)/73)
    #if (count==2) or (count==3):
    ax.set_xlabel(r'$\theta$ [deg]')
    if (count==0):# or (count==2):
        ax.set_ylabel(r'$\phi$ [deg]')
    ax.set_title(pols[count])    

plt.tight_layout()
fig.colorbar(im, ax=axes.ravel().tolist(), use_gridspec=True)
plt.savefig('amplitudes_pol.png')



fig, axes = plt.subplots(nrows=1, ncols=4, sharex = True, sharey = True, figsize = (6,4))
for count, ax in enumerate(axes.flat):
    amp_azel = np.angle(azel)
    # db = 10*np.log((amp_azel)/np.max(amp_azel)) # power dB not field
    im = ax.imshow(np.angle(polarr[:,:,count]), vmin = -np.pi, vmax = np.pi, origin = 'lower', extent = [0,180, 0, 6], aspect = (30*501)/73)
    #if (count==2) or (count==3):
    ax.set_xlabel(r'$\theta$ [deg]')
    if (count==0):# or (count==2):
        ax.set_ylabel(r'$\phi$ [deg]')
    ax.set_title(pols[count])    

plt.tight_layout()
fig.colorbar(im, ax=axes.ravel().tolist(), use_gridspec=True)
plt.savefig('phases_pol.png')




azel_types = [np.absolute(azel), np.imag(azel),np.real(azel) ,np.angle(azel)]
azel_names = ['amp', 'imag', 'real', 'phase']

for i in range(4):
    azel_type = azel_types[i]
    print (azel_type.shape)
    plt.clf()
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex = True, sharey = True, figsize = (8,6))
    for count, ax in enumerate(axes.flat):
        if i==0:
            print (i)
            azel_type[count,:,:] = 10*np.log((azel_type[count,:,:])/np.max(azel_type[count,:,:])) # power not voltage
            print (azel_type)
            print (azel_names[i])
        norm = None
        if (i==1) or (i==2):
            norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,   vmin=-1.0, vmax=1.0)
            

  
        vmax  = np.max(azel_type[0,:,:])    
        vmin  = np.min(azel_type[0,:,:])       
      
        im = ax.imshow((azel_type[count,:,:]), vmin = vmin, vmax = vmax, origin = 'lower', extent = [-6, 6, -6, 6], norm=norm)
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
    plt.savefig('azel_%s.png'%azel_names[i])


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
CDELT1  =                   %s / [deg] Coordinate increment at reference point  
CDELT2  =                   %s / [deg] Coordinate increment at reference point  
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
 """%(beamsize,beamsize,beamsize/2,(beamsize/2)+1, -(np.max(th)/beamsize), (np.max(th)/beamsize), (freq*1000)), sep='\n')

print (hdr)

















        #sin theta is radius
        # phi is angle around 
