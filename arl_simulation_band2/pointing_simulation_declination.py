"""Simulation of the effect of pointing errors on MID observations

This measures the effect of pointing errors on the change in a dirty image induced by pointing errors:
    - The pointing errors can be random per integration, static, or global, or drawn from power spectra
    - The sky can be a point source at the half power point or a realistic sky constructed from S3-SEX catalog.
    - The observation is by MID over a range of hour angles
    - Processing can be divided into chunks of time (default 1800s)
    - Dask is used to distribute the processing over a number of workers.
    - Various plots are produced, The primary output is a csv file containing information about the statistics of
    the residual images.

"""
import csv
import os
import sys
import time

from rascil.data_models.parameters import rascil_path

results_dir = rascil_path('test_results')

import numpy

from astropy.coordinates import SkyCoord
from astropy import units as u

from rascil.data_models.polarisation import PolarisationFrame

from rascil.processing_components.image.operations import show_image
from processing_library.image.operations import create_image
from rascil.processing_components.imaging.primary_beams import create_pb
from rascil.processing_components.skycomponent.operations import apply_beam_to_skycomponent
from rascil.processing_components.skycomponent.base import copy_skycomponent

import logging

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

import pprint

pp = pprint.PrettyPrinter()

# Get command line inputs
import argparse

parser = argparse.ArgumentParser(description='Simulate pointing errors')
parser.add_argument('--context', type=str, default='singlesource',
                    help='s3sky or singlesource')

parser.add_argument('--flux_limit', type=float, default=0.001, help='Flux limit (Jy)')
parser.add_argument('--show', type=str, default='False', help='Show images?')
parser.add_argument('--pbradius', type=float, default=2.0, help='Radius of sources to include (in HWHM)')
parser.add_argument('--pbtype', type=str, default='MID', help='Primary beam model: MID or MID_GAUSS')
parser.add_argument('--use_agg', type=str, default="True", help='Use Agg matplotlib backend?')
parser.add_argument('--pointing_directory', type=str, default='../../pointing_error_models/PSD_data/',
                    help='Location of pointing files')
parser.add_argument('--shared_directory', type=str, default='../../shared/',
                    help='Location of pointing files')

args = parser.parse_args()

use_agg = args.use_agg == "True"
if use_agg:
    import matplotlib as mpl
    
    mpl.use('Agg')
from matplotlib import pyplot as plt

pbtype = args.pbtype
pbradius = args.pbradius
flux_limit = args.flux_limit
shared_directory = args.shared_directory

show = args.show == 'True'
basename = os.path.basename(os.getcwd())

time_started = time.time()

# Set up details of simulated observation
nfreqwin = 1
diameter = 15.0
frequency = [1.4e9]
channel_bandwidth = [1e7]

# We need the HWHM of the primary beam. Got this by trial and error
if pbtype == 'MID':
    HWHM_deg = 0.596 * 1.4e9 / frequency[0]
elif pbtype == 'MID_GRASP':
    HWHM_deg = 0.751 * 1.4e9 / frequency[0]
elif pbtype == 'MID_GAUSS':
    HWHM_deg = 0.766 * 1.4e9 / frequency[0]
else:
    HWHM_deg = 0.596 * 1.4e9 / frequency[0]

HWHM = HWHM_deg * numpy.pi / 180.0

FOV_deg = 5.0 * HWHM_deg
print('%s: HWHM beam = %g deg' % (pbtype, HWHM_deg))

pb_npixel = 1024
d2r = numpy.pi / 180.0
pb_cellsize = d2r * FOV_deg / pb_npixel

results = list()
ra = 15.0
for declination in numpy.arange(15.0, -90.0, -15.0):
    result = dict()
    result['ra'] = ra
    result['dec'] = declination
    
    # Make a skymodel from S3
    max_flux = 0.0
    total_flux = 0.0
    print("Constructing s3sky components")
    from rascil.processing_components.simulation.testing_support import create_test_skycomponents_from_s3
    
    phasecentre = SkyCoord(ra=ra * u.deg, dec=declination * u.deg, frame='icrs', equinox='J2000')
    
    original_components = create_test_skycomponents_from_s3(flux_limit=flux_limit / 100.0,
                                                            phasecentre=phasecentre,
                                                            polarisation_frame=PolarisationFrame("stokesI"),
                                                            frequency=numpy.array(frequency),
                                                            radius=pbradius * HWHM)
    result['number_components_before'] = len(original_components)
    print("%d components before application of primary beam" % (len(original_components)))
    
    pbmodel = create_image(npixel=pb_npixel, cellsize=pb_cellsize,
                           phasecentre=phasecentre,
                           frequency=frequency,
                           polarisation_frame=PolarisationFrame("stokesI"))
    # Use MID_GAUSS to filter the components since MID_GRASP is in local coordinates
    pb = create_pb(pbmodel, "MID_GAUSS", pointingcentre=phasecentre, use_local=False)
    pb_applied_components = [copy_skycomponent(c) for c in original_components]
    pb_applied_components = apply_beam_to_skycomponent(pb_applied_components, pb)
    filtered_components = []
    for icomp, comp in enumerate(pb_applied_components):
        if comp.flux[0, 0] > flux_limit:
            total_flux += comp.flux[0, 0]
            if abs(comp.flux[0, 0]) > max_flux:
                max_flux = abs(comp.flux[0, 0])
            filtered_components.append(original_components[icomp])
    result['flux_limit'] = flux_limit
    result['number_components'] = len(filtered_components)
    result['max_flux'] = max_flux
    result['total_flux'] = total_flux
    print("RA %.1f Dec %.1f number components %d peak flux %.3f total flux %.3f" % (ra, declination,
                                                                                    len(filtered_components),
                                                                                    max_flux, total_flux))
    original_components = [copy_skycomponent(c) for c in filtered_components]
    plt.clf()
    show_image(pb, components=original_components)
    plt.savefig('component_image_flux_limit_%.3f_dec_%.1f.png' % (flux_limit, declination))
    plt.show(block=False)
    results.append(result)

plt.clf()
mflux = [r['max_flux'] for r in results]
tflux = [r['total_flux'] for r in results]
dec = [r['dec'] for r in results]
plt.plot(dec, mflux, '.', label='Max flux')
plt.plot(dec, tflux, '.', label='Total flux')
plt.xlabel('Declination (deg)')
plt.ylabel('Flux')
plt.legend()
plt.savefig('fluxes_vs_declination_%.3fJy.png' % (flux_limit))
plt.show(block=False)

plt.clf()
ncomps = [r['number_components'] for r in results]
nbcomps = [r['number_components_before'] for r in results]
dec = [r['dec'] for r in results]
plt.semilogy(dec, nbcomps, '.', label='Number before PB')
plt.semilogy(dec, ncomps, '.', label='Number after PB')
plt.xlabel('Declination (deg)')
plt.ylabel('Number components')
plt.legend()
plt.savefig('number_components_vs_declination_%.3fJy.png' % (flux_limit))
plt.show(block=False)

pp.pprint(results)

filename = 'model_sky_%.1f.csv' % flux_limit
with open(filename, 'a') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=results[0].keys(), delimiter=',', quotechar='|',
                            quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    for result in results:
        writer.writerow(result)
    csvfile.close()
