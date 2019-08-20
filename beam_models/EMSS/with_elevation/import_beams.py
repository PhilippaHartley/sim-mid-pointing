#!/usr/bin/env python3

import matplotlib as mpl
import numpy as np
import scipy.interpolate
import numpy.fft
import scipy.io

mpl.use('agg')
import matplotlib.pyplot as plt

import numpy

from astropy.wcs import WCS

from data_models.polarisation import PolarisationFrame

from processing_library.image.operations import create_image_from_array
from processing_components.image.operations import export_image_to_fits


def create_arl_image(pol_planes, cellsize, frequency, channel_bandwidth=1e6, shift_peak=True):
    ny, nx = pol_planes[0].shape
    
    assert len(pol_planes) == 4
    
    w = WCS(naxis=4)
    # The negation in the longitude is needed by definition of RA, DEC
    w.wcs.cdelt = [-cellsize, cellsize, -1.0, channel_bandwidth]
    w.wcs.crpix = [nx // 2, ny // 2, 1.0, 1.0]
    w.wcs.ctype = ['AZELGEO long', 'AZELGEO lati', 'STOKES', 'FREQ']
    w.wcs.crval = [0.0, 0.0, -5.0, frequency]
    w.wcs.cunit = ['deg', 'deg', '', 'Hz']
    w.naxis = 4
    w.wcs.radesys = 'ICRS'
    w.wcs.equinox = 2000.0
    
    beam_out = numpy.zeros([1, 4, ny, nx], dtype=complex)
    for pol in range(4):
        beam_out[0, pol, ...] = numpy.transpose(pol_planes[pol])
    
    # The following transforms are guessed to provide realistic looking beams
    # 1. Renormalise
    beam_out /= numpy.max(numpy.abs(beam_out))
    # 2. Remove phase error in image plane
    beam_out *= numpy.conjugate(beam_out[0, 0, ny // 2, nx // 2])

    # 3. Remove phase gradient in image plane
    dy = numpy.mod(numpy.angle(beam_out[0, 0, ny//2 + 1, nx // 2]) -
                   numpy.angle(beam_out[0, 0, ny//2 - 1, nx // 2]), numpy.pi) / 2.0
    rotator = numpy.exp(-1.0j * dy * (numpy.arange(ny) - ny / 2.0))
    beam_out *= rotator[numpy.newaxis, numpy.newaxis, :, numpy.newaxis]
    for pol in [1, 3]:
        beam_out[:, pol, ...] = -1.0 * beam_out[:, pol, ...]
        
    # FFT interpolate to get peak on nx//2, ny//2
    if shift_peak:
        beam_xfr = numpy.fft.fft2(beam_out, axes=(0,1))
        print(numpy.angle(beam_xfr[0, 0, ny//2 + 1, nx // 2]), numpy.angle(beam_xfr[0, 0, ny//2 - 1, nx // 2]))
        dv = (numpy.angle(beam_xfr[0, 0, ny//2 + 1, nx // 2]) - numpy.angle(beam_xfr[0, 0, ny//2 - 1, nx // 2])) / 2.0
        print(dv)
        rotator = numpy.exp(1.0j * dv * (numpy.arange(ny) - ny // 2))
        beam_xfr *= rotator[numpy.newaxis, numpy.newaxis, :, numpy.newaxis]
        beam_out = numpy.fft.ifft2(beam_xfr, axes=(0,1))

    vp_real = create_image_from_array(beam_out.real, w, polarisation_frame=PolarisationFrame("linear"))
    vp_imag = create_image_from_array(beam_out.imag, w, polarisation_frame=PolarisationFrame("linear"))
    vp_amp = create_image_from_array(numpy.abs(beam_out), w, polarisation_frame=PolarisationFrame("linear"))
    vp_phase = create_image_from_array(numpy.angle(beam_out), w, polarisation_frame=PolarisationFrame("linear"))
    
    return vp_real, vp_imag, vp_amp, vp_phase


def interpolate_beam(th, ph, beam_inp, n, extent):
    nx, ny = n
    xmin, xmax, ymin, ymax = extent
    
    # Set x and y values of output grid.
    
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    
    # Make x and y into 2D grids and flatten into 1D arrays.
    
    x2, y2 = np.meshgrid(x, y)
    x2.shape = (-1,)
    y2.shape = (-1,)
    
    # Get (th, ph) values of (x, y) grid points assuming
    #
    # x = th * cos(ph)
    # y = th * sin(ph)
    
    thi = np.sqrt(x2 ** 2 + y2 ** 2)
    phi = np.mod(np.degrees(np.arctan2(y2, x2)), 360.0)
    xi = np.stack((thi, phi), axis=1)
    
    # Interpolate real and imaginary parts separately then combine.
    
    tmp_real = scipy.interpolate.interpn((th, ph), beam_inp.real, xi)
    tmp_imag = scipy.interpolate.interpn((th, ph), beam_inp.imag, xi)
    beam_out = tmp_real + 1j * tmp_imag
    
    # Reshape output into 2D image.
    
    beam_out.shape = (nx, ny)
    
    return beam_out


def plot_beam(beam, title, extent):
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.suptitle(title)
    
    ax = axes[0, 0]
    ax.set_title('real')
    ax.set_ylabel('y / deg')
    im = ax.imshow(beam.real, extent=extent, origin='lower')
    fig.colorbar(im, ax=ax)
    
    ax = axes[0, 1]
    ax.set_title('imag')
    im = ax.imshow(beam.imag, extent=extent, origin='lower')
    fig.colorbar(im, ax=ax)
    
    ax = axes[1, 0]
    ax.set_xlabel('x / deg')
    ax.set_ylabel('y / deg')
    ax.set_title('amplitude')
    im = ax.imshow(np.abs(beam), extent=extent, origin='lower')
    fig.colorbar(im, ax=ax)
    
    ax = axes[1, 1]
    ax.set_title('phase')
    ax.set_xlabel('x / deg')
    im = ax.imshow(np.angle(beam, deg=True), extent=extent, origin='lower')
    fig.colorbar(im, ax=ax)


iform = '{b}_{e}_{f}.mat'
oform = '{b}_{e}_{f:04d}_{j}.png'
tform = 'band = {b}, elev = {e} deg, freq = {f} MHz, jones = {j}'
fitsform = '{b}_{e}_{f:04d}_{t}_shift.fits'

# n = nx, ny
# extent = xmin, xmax, ymin, ymax

n = 1001, 1001
n = 281, 281
extent = -4.0, 4.0, -4.0, 4.0
cellsize = 8 / n[0]

# band = [('B1', [365, 415, 465, 515, 565, 615, 665, 715, 765, 815, 865,
#                915, 965, 1015, 1050]),
#        ('B2', [965, 1000, 1060, 110, 1160,
#                1220, 1252, 1310, 1360, 1410, 1460, 1510, 1610, 1660, 1710,
#                1760]),
#        ('Ku', [11452, 11697, 11699, 11700, 12179, 12251, '12501_5'])]

band = [('B2', [965, 1000, 1060, 1100, 1160, 1220, 1252, 1310, 1360,
                1410, 1460, 1510, 1610, 1660, 1710, 1760])]
band = [('B2', [1360])]

elev = [15, 45, 90]
jones = ['Jpv', 'Jqh', 'Jph', 'Jqv']

for b, freq in band:
    for e in elev:
        for f in freq:
            
            print(b, e, f)
            ifile = iform.format(b=b, e=e, f=f)
            data = scipy.io.loadmat(ifile)
            ph = data['ph'].squeeze()
            th = data['th'].squeeze()
            
            pol_planes = list()
            for j in jones:
                print(b, e, f, j)
                ofile = oform.format(b=b, e=e, f=f, j=j)
                title = tform.format(b=b, e=e, f=f, j=j)
                beam_inp = data[j]
                beam_out = interpolate_beam(th, ph, beam_inp, n, extent)
                pol_planes.append(beam_out)
                plot_beam(beam_out, title, extent)
                plt.savefig(ofile)
                plt.close()
            
            vp_real, vp_imag, vp_amp, vp_phase = create_arl_image(pol_planes, cellsize, f, shift_peak=True)
            export_image_to_fits(vp_real, fitsform.format(b=b, e=e, f=f, t='real'))
            export_image_to_fits(vp_imag, fitsform.format(b=b, e=e, f=f, t='imag'))
            export_image_to_fits(vp_amp, fitsform.format(b=b, e=e, f=f, t='amp'))
            export_image_to_fits(vp_phase, fitsform.format(b=b, e=e, f=f, t='phase'))
