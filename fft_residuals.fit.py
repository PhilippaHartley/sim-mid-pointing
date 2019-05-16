
import os
import sys

import numpy

from matplotlib import pylab

pylab.rcParams['figure.figsize'] = (8.0, 8.0)
pylab.rcParams['image.cmap'] = 'rainbow'

from matplotlib import pyplot as plt
from processing_library.image.operations import fft_image
from wrappers.serial.image.operations import show_image, import_image_from_fits, export_image_to_fits, \
    qa_image

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

for sim in range(4):

    pes = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]

    for pe in pes:
        residual_image = import_image_from_fits('simulation%d/PE_%.1f_arcsec_arl.fits' % (sim, pe))
        residual_fft = fft_image(residual_image)
        assert numpy.max(numpy.abs(numpy.imag(residual_fft.data))) < 1e-15
        
        residual_fft.data = numpy.real(residual_fft.data)
        print(qa_image(residual_fft, context='simulation%d/PE_%.1f_arcsec_arl.fits' % (sim, pe)))
        show_image(residual_fft, title='simulation%d/PE_%.1f_arcsec_arl_fft.fits' % (sim, pe))
        plt.savefig('simulation%d/PE_%.1f_arcsec_arl_fft.png' % (sim, pe) )
        plt.show()
        export_image_to_fits(residual_fft, 'simulation%d/PE_%.1f_arcsec_arl_fft.fits' % (sim, pe))
