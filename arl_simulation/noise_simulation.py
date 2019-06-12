"""Simulation of the effect of pointing errors on MID observations
"""
import csv
import socket
import sys
import time

import seqfile

from data_models.parameters import arl_path

results_dir = arl_path('test_results')

import numpy

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

import matplotlib as mpl
mpl.use('Agg')

from matplotlib import pyplot as plt

from data_models.polarisation import PolarisationFrame
from data_models.memory_data_models import Skycomponent, SkyModel
from data_models.data_model_helpers import export_pointingtable_to_hdf5

from wrappers.serial.visibility.base import create_blockvisibility
from wrappers.serial.image.operations import show_image, qa_image, export_image_to_fits
from wrappers.serial.simulation.configurations import create_configuration_from_MIDfile
from wrappers.serial.simulation.testing_support import simulate_pointingtable
from wrappers.serial.imaging.primary_beams import create_vp, create_pb
from wrappers.serial.imaging.base import create_image_from_visibility, advise_wide_field
from wrappers.serial.calibration.pointing import create_pointingtable_from_blockvisibility
from wrappers.serial.simulation.pointing import create_gaintable_from_pointingtable
from wrappers.arlexecute.visibility.base import copy_visibility
from wrappers.arlexecute.visibility.coalesce import convert_blockvisibility_to_visibility, \
    convert_visibility_to_blockvisibility


from workflows.arlexecute.skymodel.skymodel_arlexecute import predict_skymodel_list_compsonly_arlexecute_workflow
from workflows.arlexecute.imaging.imaging_arlexecute import invert_list_arlexecute_workflow
from workflows.serial.imaging.imaging_serial import weight_list_serial_workflow

from wrappers.arlexecute.execution_support.arlexecute import arlexecute
from wrappers.arlexecute.execution_support.dask_init import get_dask_Client

from wrappers.serial.image.operations import import_image_from_fits
from wrappers.serial.simulation.ionospheric_screen import create_gaintable_from_screen

from processing_components.simulation.noise import create_gaintable_from_noise_sources

import logging

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

if __name__ == '__main__':
    
    # Get command line inputs
    import argparse
    
    parser = argparse.ArgumentParser(description='Simulate pointing errors')
    parser.add_argument('--context', type=str, default='singlesource',
                        help='s3sky or singlesource')
    
    parser.add_argument('--rmax', type=float, default=1e5,
                        help='Maximum distance of station from centre (m)')
    
    parser.add_argument('--global_pe', type=float, nargs=2, default=[0.0, 0.0], help='Global pointing error')
    parser.add_argument('--static_pe', type=float, default=0.0, help='Multiplier for static errors')
    parser.add_argument('--dynamic_pe', type=float, default=1.0, help='Multiplier for dynamic errors')
    parser.add_argument('--nnodes', type=int, default=1, help='Number of nodes')
    parser.add_argument('--nthreads', type=int, default=1, help='Number of threads')
    parser.add_argument('--memory', type=int, default=8, help='Memory per worker')
    parser.add_argument('--nworkers', type=int, default=8, help='Number of workers')
    parser.add_argument('--flux_limit', type=float, default=1.0, help='Flux limit (Jy)')
    parser.add_argument('--show', type=str, default='False', help='Show images?')
    parser.add_argument('--ngroup', type=int, default=8, help='Process in groups this large')
    parser.add_argument('--outlier', type=str, default='False', help='Do the outlier field as well?')
    parser.add_argument('--seed', type=int, default=18051955, help='Random number seed')
    parser.add_argument('--snapshot', type=str, default=False, help='Do snapshot only')

    args = parser.parse_args()
    
    snapshot = args.snapshot == 'True'
    
    seed = args.seed
    
    print("Random number seed is ", seed)
    show = args.show == 'True'
    outlier = args.outlier == 'True'
    context = args.context
    rmax = args.rmax
    flux_limit = args.flux_limit
    
    nworkers = args.nworkers
    nnodes = args.nnodes
    threads_per_worker = args.nthreads
    memory = args.memory
    
    ngroup = args.ngroup
    
    print("Using %s workers" % nworkers)
    print("Using %s threads per worker" % threads_per_worker)
    
    # Set up DASK
    client = get_dask_Client(threads_per_worker=threads_per_worker,
                             processes=threads_per_worker == 1,
                             memory_limit=memory * 1024 * 1024 * 1024,
                             n_workers=nworkers)
    arlexecute.set_client(client=client)
    
    # Set up details of simulated observation
    nfreqwin = 1
    ntimes = 10
    diameter = 15.0
    frequency = [1.4e9]
    channel_bandwidth = [1e7]
    
    if snapshot:
        ntimes = 1
        times = [0.0]
    else:
        
        h2r = numpy.pi / 12.0
        times = numpy.linspace(-6 * h2r, +6 * h2r, ntimes)
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    outlier_phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-35.0 * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    mid = create_configuration_from_MIDfile('../shared/ska1mid_small.cfg', rmax=rmax, location=location)

    block_vis = create_blockvisibility(mid, times, frequency=frequency,
                                           channel_bandwidth=channel_bandwidth, weight=1.0,
                                           phasecentre=phasecentre,
                                           polarisation_frame=PolarisationFrame("stokesI"), zerow=False)
    print(block_vis)
    vis = convert_blockvisibility_to_visibility(block_vis)
    
    vis.data['uvw'][...,2] = 0.0
    
    # We need the HWHM of the primary beam. Got this by trial and error
    HWHM_deg = 1.03 * 180.0 * 3e8 / (numpy.pi * diameter * frequency[0])
    
    print('HWHM beam = %g deg' % HWHM_deg)
    HWHM = HWHM_deg * numpy.pi / 180.0

    advice = advise_wide_field(vis, guard_band_image=1.0, delA=0.02)

    cellsize = advice['cellsize']
    if context == 's3sky':
        pb_npixel = 4096
        pb_cellsize = 4.0 * HWHM / pb_npixel
        npixel = pb_npixel
    else:
        npixel = 512
        pb_npixel = 4096
        pb_cellsize = HWHM / pb_npixel
    
    if show:
        plt.clf()
        plt.plot(-vis.u, -vis.v, '.', color='b', markersize=0.2)
        plt.plot(vis.u, vis.v, '.', color='b', markersize=0.2)
        plt.savefig('uvcoverage.png')
        plt.show(block=False)

    # Uniform weighting
    model = create_image_from_visibility(vis, npixel=npixel, frequency=frequency,
                                         nchan=nfreqwin, cellsize=cellsize, phasecentre=phasecentre,
                                         polarisation_frame=PolarisationFrame("stokesI"))
    print(model)
    vis = weight_list_serial_workflow([vis], [model])[0]
    block_vis = convert_visibility_to_blockvisibility(vis)

    print("Inverting to get on-source PSF")
    psf = invert_list_arlexecute_workflow([vis], [model], '2d', dopsf=True)
    psf, sumwt = arlexecute.compute(psf, sync=True)[0]

    export_image_to_fits(psf, 'PSF_arl.fits')

    if show:
        show_image(psf, cm='gray_r', title='PSF', vmin=-0.01, vmax=0.1)
        plt.savefig('PSF_arl.png')
        plt.show(block=False)
        
    # Construct the skycomponents
    if context == 'singlesource':
        print("Constructing single component")
        # Put a single point source at the phasecentre
        original_components = [Skycomponent(flux=[[1.0]], direction=phasecentre,
                                            frequency=frequency,
                                            polarisation_frame=PolarisationFrame('stokesI'))]
        
        offset = [180.0 * pb_cellsize * pb_npixel / (2.0 * numpy.pi), 0.0]
        HWHM = HWHM_deg * numpy.pi / 180.0
        # The primary beam is offset to approximately the halfpower point
        pb_direction = SkyCoord(ra=(+15.0 + offset[0] / numpy.cos(-45.0 * numpy.pi / 180.0)) * u.deg,
                                dec=(-45.0 + offset[1]) * u.deg, frame='icrs', equinox='J2000')
    
    else:
        # Make a skymodel from S3
        print("Constructing s3sky components")
        from wrappers.serial.simulation.testing_support import create_test_skycomponents_from_s3
        
        original_components = create_test_skycomponents_from_s3(flux_limit=flux_limit,
                                                                phasecentre=phasecentre,
                                                                polarisation_frame=PolarisationFrame("stokesI"),
                                                                frequency=numpy.array(frequency),
                                                                radius=pb_cellsize * pb_npixel)
        HWHM = HWHM_deg * numpy.pi / 180.0
        # Primary beam points to the phasecentre
        pb_direction = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    
    # ### Calculate the voltage patterns with and without pointing errors
    vp = create_image_from_visibility(block_vis, npixel=pb_npixel, frequency=frequency,
                                      nchan=nfreqwin, cellsize=pb_cellsize, phasecentre=phasecentre,
                                      override_cellsize=False)
    
    pb = create_pb(vp, 'MID', pointingcentre=pb_direction)
    
    if show:
        show_image(pb, title='%s: primary beam' % context)
        plt.savefig('PB_arl.png')
        plt.show(block=False)
    
    print("Constructing voltage pattern")
    vp = create_vp(vp, 'MID', pointingcentre=pb_direction)
    pt = create_pointingtable_from_blockvisibility(block_vis, vp)
    
    no_error_pt = simulate_pointingtable(pt, 0.0, 0.0, seed=seed)
   # export_pointingtable_to_hdf5(no_error_pt, 'pointingsim_%s_noerror_pointingtable.hdf5' % context)
    error_noise = create_gaintable_from_noise_sources(block_vis, original_components, no_error_pt, vp)


    
    no_error_sm = [SkyModel(components=[original_components[i]], gaintable=error_noise[i])
                   for i, _ in enumerate(original_components)]
    
    # We do this in chunks of eight to avoid creating all visibilities at once
    no_error_blockvis = copy_visibility(block_vis, zero=True)
    
    print("Predicting error-free visibilities in chunks of %d skymodels" % ngroup)
    future_vis = arlexecute.scatter(no_error_blockvis)
    chunks = [no_error_sm[i:i + ngroup] for i in range(0, len(no_error_sm), ngroup)]
    for chunk in chunks:
        temp_vis = predict_skymodel_list_compsonly_arlexecute_workflow(future_vis, chunk, context='2d', docal=True)
        work_vis = arlexecute.compute(temp_vis, sync=True)
        for w in work_vis:
            no_error_blockvis.data['vis'] += w.data['vis']
        assert numpy.max(numpy.abs(no_error_blockvis.data['vis'])) > 0.0
    
    no_error_vis = convert_blockvisibility_to_visibility(no_error_blockvis)

    dirty = invert_list_arlexecute_workflow([no_error_vis], [model], '2d')
    dirty, sumwt = arlexecute.compute(dirty, sync=True)[0]

    export_image_to_fits(dirty, 'with_noise.fits')
    
    pes = [1.0]#, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
    results = []
    
    filename = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.csv')
    print('Saving results to %s' % filename)
    plotfile = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.jpg')
    
    epoch = time.strftime("%Y-%m-%d %H:%M:%S")
    
