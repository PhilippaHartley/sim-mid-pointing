
# coding: utf-8

# # Imaging and deconvolution demonstration

# This script makes a fake data set and then deconvolves it. Finally the full and residual visibility are plotted.

# In[1]:


import os
import sys

sys.path.append(os.path.join('..', '..'))

from rascil.data_models.parameters import rascil_path
results_dir = rascil_path('test_results')


from matplotlib import pylab

pylab.rcParams['figure.figsize'] = (8.0, 8.0)
pylab.rcParams['image.cmap'] = 'rainbow'

import numpy

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord

from matplotlib import pyplot as plt

from processing_library.util.coordinate_support import *

from rascil.processing_components.image.iterators import image_raster_iter

from rascil.processing_components.visibility.base import create_visibility
from rascil.processing_components.skycomponent.operations import create_skycomponent
from rascil.processing_components.image.operations import show_image, export_image_to_fits
from rascil.processing_components.image.deconvolution import deconvolve_cube, restore_cube
from rascil.processing_components.visibility.iterators import vis_timeslice_iter
from rascil.processing_components.simulation.configurations import create_named_configuration
from rascil.processing_components.simulation.testing_support import create_test_image
from rascil.processing_components.imaging.base import create_image_from_visibility
from rascil.processing_components.imaging.base import advise_wide_field

from rascil.workflows.serial.imaging.imaging_serial import invert_list_serial_workflow, predict_list_serial_workflow

from rascil.data_models.polarisation import PolarisationFrame

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

mpl_logger = logging.getLogger("matplotlib") 
mpl_logger.setLevel(logging.WARNING) 


# In[2]:


pylab.rcParams['figure.figsize'] = (6.0, 6.0)
pylab.rcParams['image.cmap'] = 'rainbow'


# Construct LOW core configuration

# In[3]:


def plotuvw(uvwcoord, x=0, y=1, z=2):
    plt.clf()
    plt.plot(uvwcoord[:,x], uvwcoord[:,y], '.', color='b')
    plt.show()


# In[5]:


def xyz_to_uvw(xyz, ha, dec):
    """
    Rotate :math:`(x,y,z)` positions in earth coordinates to
    :math:`(u,v,w)` coordinates relative to astronomical source
    position :math:`(ha, dec)`. Can be used for both antenna positions
    as well as for baselines.

    Hour angle and declination can be given as single values or arrays
    of the same length. Angles can be given as radians or astropy
    quantities with a valid conversion.

    :param xyz: :math:`(x,y,z)` co-ordinates of antennas in array
    :param ha: hour angle of phase tracking centre (:math:`ha = ra - lst`)
    :param dec: declination of phase tracking centre.
    """
    
    x, y, z = numpy.hsplit(xyz, 3)  # pylint: disable=unbalanced-tuple-unpacking
    
    # Two rotations:
    #  1. by 'ha' along the z axis
    #  2. by '90-dec' along the u axis
    u = x * numpy.cos(ha) - y * numpy.sin(ha)
    v0 = x * numpy.sin(ha) + y * numpy.cos(ha)
    w = z * numpy.sin(dec) - v0 * numpy.cos(dec)
    lat = 180.0 * numpy.arctan2(v0,z) / numpy.pi
    v = z * numpy.cos(dec) + v0 * numpy.sin(dec)
    
    return numpy.hstack([u, v, w])


# In[7]:


antfile=rascil_path("data/configurations/MID_SKA-TEL-INSA-0000537_Rev05.txt")
antdiamlonglat = numpy.genfromtxt(antfile, usecols=[0, 1, 2], delimiter="\t")
location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)

print(location.geocentric)
    
assert antdiamlonglat.shape[1] == 3, ("Antenna array has wrong shape %s" % antdiamlonglat.shape)
xyz = numpy.zeros([antdiamlonglat.shape[0] - 1, 3])
for ant in range(antdiamlonglat.shape[0] - 1):
    loc = EarthLocation(lon=antdiamlonglat[ant, 1], lat=antdiamlonglat[ant, 2], height=0.0).geocentric
    xyz[ant] = [loc[0].to(u.m).value, loc[1].to(u.m).value, loc[2].to(u.m).value]
    
uvw=xyz_to_uvw(xyz, -location.geodetic[0].to('rad').value, location.geodetic[1].to('rad').value)

print(uvw[::19,:])

