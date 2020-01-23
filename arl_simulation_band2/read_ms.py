import os
import sys

sys.path.append(os.path.join('..', '..'))

from rascil.data_models.parameters import rascil_path

results_dir = rascil_path('test_results')

from matplotlib import pylab

pylab.rcParams['figure.figsize'] = (8.0, 8.0)
pylab.rcParams['image.cmap'] = 'rainbow'

import numpy

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants

from matplotlib import pyplot as plt

from rascil.processing_components.visibility.base import create_visibility
from rascil.processing_components.simulation.configurations import create_named_configuration

from rascil.data_models.polarisation import PolarisationFrame

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

pylab.rcParams['figure.figsize'] = (12.0, 12.0)
pylab.rcParams['image.cmap'] = 'rainbow'

from rascil.processing_components.visibility.base import create_blockvisibility_from_ms

vis = create_blockvisibility_from_ms('PE_0.0_arcsec.ms')[0]

print(vis.uvw.shape)
plt.clf()
uvw=vis.uvw.reshape([65*197*197,3])
vu=uvw[...,0]
vv=uvw[...,1]
plt.plot(vu, vv, '.', color='b')
plt.plot(-vu, -vv, '.', color='b')
plt.show()

print(vis.configuration)