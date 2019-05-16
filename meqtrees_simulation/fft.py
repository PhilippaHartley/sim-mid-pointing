import scipy

import numpy as np

from scipy.fftpack import fft2

from astropy.io import fits
from matplotlib import pyplot as plt


tb.open('simulation30/setupvis.ms')
res = tb.getcol('MODEL_DATA')
uvw = tb.getcol('UVW')
print res[0][0].shape
print uvw.shape
print uvw[0].shape
np

plt.hist(np.real(res[0][0]), bins = 50, alpha = 0.5)
plt.hist(np.imag(res[0][0]), bins = 50, alpha = 0.5)
plt.savefig('res_hist.png')
plt.clf()
plt.scatter(np.hypot(uvw[0],uvw[1]), np.real(res[0][0]))




arr = np.zeros((    np.max(uvw[0])+1  ,   np.max(uvw[1])  +1 ))
print arr.shape
for i in range(len(res[0][0])):
    if np.int(uvw[0][i]) > 0:
        if np.int(uvw[1][i]) >0:

            arr[np.int(uvw[0][i]),np.int(uvw[1][i])]+=np.real(res[0][0][i])

fits.writeto('mssetupreal.fits', arr[:1000,:1000])

print arr

plt.clf()
plt.imshow(arr[:1000,:1000])
plt.savefig('arr.png')





plt.savefig('uvw.png')
plt.clf()
a = fits.getdata('simulation30/residuals_1_arcsec_normal.fits')
a = a.astype(np.float)
plt.imshow(a)
plt.savefig('res.png')
b = fft2(a)
#fits.writeto('real.fits', np.real(b))
#f#its.writeto('imag.fits', np.imag(b))
plt.imshow(np.real(b))
plt.savefig('real.png')
plt.imshow(np.imag(b))
plt.savefig('imag.png')

