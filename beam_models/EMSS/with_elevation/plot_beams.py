#!/usr/bin/env python3

import numpy as np
import scipy.io
import scipy.interpolate
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt


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

    thi = np.sqrt(x2**2 + y2**2)
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

# n = nx, ny
# extent = xmin, xmax, ymin, ymax

n = 1001, 1001
extent = -4.0, 4.0, -4.0, 4.0

#band = [('B1', [365, 415, 465, 515, 565, 615, 665, 715, 765, 815, 865,
#                915, 965, 1015, 1050]),
#        ('B2', [965, 1000, 1060, 110, 1160,
#                1220, 1252, 1310, 1360, 1410, 1460, 1510, 1610, 1660, 1710,
#                1760]),
#        ('Ku', [11452, 11697, 11699, 11700, 12179, 12251, '12501_5'])]

band = [('B2', [965, 1000, 1060, 1100, 1160, 1220, 1252, 1310, 1360,
                1410, 1460, 1510, 1610, 1660, 1710, 1760])]

elev = [15, 45, 90]
jones = ['Jph', 'Jpv', 'Jqh', 'Jqv']
                
for b, freq in band:
    for e in elev:
        for f in freq:
                
            print(b, e, f)
            ifile = iform.format(b=b, e=e, f=f)
            data = scipy.io.loadmat(ifile)
            ph = data['ph'].squeeze()
            th = data['th'].squeeze()
                
            for j in jones:
                
                print(b, e, f, j)
                ofile = oform.format(b=b, e=e, f=f, j=j)
                title = tform.format(b=b, e=e, f=f, j=j)
                beam_inp = data[j]
                beam_out = interpolate_beam(th, ph, beam_inp, n, extent)
                plot_beam(beam_out, title, extent)
                plt.savefig(ofile)
                plt.close()
