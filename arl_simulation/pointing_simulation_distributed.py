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
import socket
import sys
import time

import seqfile

from data_models.parameters import arl_path

results_dir = arl_path('test_results')

import numpy

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from data_models.polarisation import PolarisationFrame
from data_models.memory_data_models import Skycomponent, SkyModel

from processing_library.image.operations import create_empty_image_like
from wrappers.serial.visibility.base import create_blockvisibility
from wrappers.serial.image.operations import show_image, qa_image, export_image_to_fits
from wrappers.serial.simulation.configurations import create_configuration_from_MIDfile
from wrappers.serial.simulation.testing_support import simulate_pointingtable, simulate_pointingtable_from_timeseries
from wrappers.serial.imaging.primary_beams import create_vp, create_pb
from wrappers.serial.imaging.base import create_image_from_visibility, advise_wide_field
from wrappers.serial.calibration.pointing import create_pointingtable_from_blockvisibility
from wrappers.serial.simulation.pointing import simulate_gaintable_from_pointingtable
from processing_library.util.coordinate_support import hadec_to_azel
from wrappers.arlexecute.visibility.base import copy_visibility
from wrappers.arlexecute.visibility.coalesce import convert_blockvisibility_to_visibility, \
    convert_visibility_to_blockvisibility

from workflows.arlexecute.skymodel.skymodel_arlexecute import predict_skymodel_list_compsonly_arlexecute_workflow
from workflows.arlexecute.imaging.imaging_arlexecute import invert_list_arlexecute_workflow
from workflows.arlexecute.imaging.imaging_arlexecute import weight_list_arlexecute_workflow
from workflows.shared.imaging.imaging_shared import sum_invert_results

from wrappers.arlexecute.execution_support.arlexecute import arlexecute
from wrappers.arlexecute.execution_support.dask_init import get_dask_Client

import logging

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

import pprint

pp = pprint.PrettyPrinter()


# Process a set of BlockVisibility's, creating pointing errors, converting to gainables, applying
# the gaintables to the FT of the skycomponents, and dirty images, one per BlockVisibility
def create_vis_list_with_errors(bvis_list, original_components, model_list, vp_list,
                                use_radec=False, pointing_error=0.0, static_pointing_error=0.0,
                                global_pointing_error=None, time_series='', seeds=None,
                                reference_pointing=False, pointing_file=None):
    if global_pointing_error is None:
        global_pointing_error = [0.0, 0.0]
    
    # One pointing table per visibility
    
    pt_list = [arlexecute.execute(create_pointingtable_from_blockvisibility)(bvis) for bvis in bvis_list]
    if time_series is '':
        pt_list = [arlexecute.execute(simulate_pointingtable)(pt, pointing_error=pointing_error,
                                                              static_pointing_error=static_pointing_error,
                                                              global_pointing_error=global_pointing_error,
                                                              seed=seeds[ipt])
                   for ipt, pt in enumerate(pt_list)]
    else:
        pt_list = [arlexecute.execute(simulate_pointingtable_from_timeseries)(pt, type=time_series,
                                                                              pointing_file=pointing_file,
                                                                              reference_pointing=reference_pointing)
                   for pt in pt_list]
    
    if show:
        tmp_error_pt_list = arlexecute.compute(pt_list, sync=True)
        plt.clf()
        r2a = 180.0 * 3600.0 / numpy.pi
        for pt in tmp_error_pt_list:
            plt.plot(pt.time, r2a * pt.pointing[:, 0, 0, 0, 0], '.', color='r')
            plt.plot(pt.time, r2a * pt.pointing[:, 0, 0, 0, 1], '.', color='b')
        plt.title("%s: dish 0 pointing" % basename)
        plt.xlabel('Time (s)')
        plt.ylabel('Offset (arcsec)')
        plt.savefig('pointing_error.png')
        plt.show(block=False)
    
    # Create the gain tables, one per Visibility and per component
    gt_list = [arlexecute.execute(simulate_gaintable_from_pointingtable)
               (bvis, original_components, pt_list[ibv], vp_list[ibv], use_radec=use_radec)
               for ibv, bvis in enumerate(bvis_list)]
    if show:
        tmp_gt_list = arlexecute.compute(gt_list, sync=True)
        plt.clf()
        for gt in tmp_gt_list:
            amp = numpy.abs(gt[0].gain[:, 0, 0, 0, 0])
            plt.semilogy(gt[0].time[amp > 0.0], 1.0 / amp[amp > 0.0], '.')
        plt.ylim([1e-3, 1.1])
        plt.title("%s: dish 0 amplitude gain" % basename)
        plt.xlabel('Time (s)')
        plt.savefig('gaintable.png')
        plt.show(block=False)
    
    # Each component in original components becomes a separate skymodel
    # Inner nest is over skymodels, outer is over bvis's
    error_sm_list = [[
        arlexecute.execute(SkyModel, nout=1)(components=[original_components[i]], gaintable=gt_list[ibv][i])
        for i, _ in enumerate(original_components)] for ibv, bv in enumerate(bvis_list)]
    
    # Predict each visibility for each skymodel. We keep all the visibilities separate
    # and add up dirty images at the end of processing. We calibrate which applies the voltage pattern
    error_bvis_list = [arlexecute.execute(copy_visibility, nout=1)(bvis, zero=True) for bvis in bvis_list]
    error_bvis_list = [predict_skymodel_list_compsonly_arlexecute_workflow(error_bvis_list[ibv], error_sm_list[ibv],
                                                                           context='2d', docal=True)
                       for ibv, bvis in enumerate(error_bvis_list)]
    
    # Inner nest is bvis per skymodels, outer is over vis's. First we convert all blockvisibilities to
    # visibilities.
    error_vis_list = [[arlexecute.execute(convert_blockvisibility_to_visibility)(bvis[i])
                       for i, _ in enumerate(original_components)] for bvis in error_bvis_list]
    
    # Now for each visibility/component, we make the dirty images
    dirty_list = list()
    # We just add the component dirty images since the weights should be the same
    def sum_images(images):
        sum_image = create_empty_image_like(images[0][0])
        for im in images:
            sum_image.data += im[0].data
        return sum_image, images[0][1]
    
    dirty_list = list()
    for vis in error_vis_list:
        assert len(vis) == len(original_components)
        result = invert_list_arlexecute_workflow(vis, model_list, '2d')
        dirty_list.append(arlexecute.execute(sum_images)(result))

    return dirty_list


if __name__ == '__main__':
    
    print(" ")
    print("Distributed simulation of pointing errors for SKA-MID")
    print("-----------------------------------------------------")
    print(" ")

    memory_use = dict()
    
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
    parser.add_argument('--memory', type=int, default=8, help='Memory per worker (GB)')
    parser.add_argument('--nworkers', type=int, default=8, help='Number of workers')
    parser.add_argument('--flux_limit', type=float, default=1.0, help='Flux limit (Jy)')
    parser.add_argument('--show', type=str, default='False', help='Show images?')
    parser.add_argument('--ngroup', type=int, default=8, help='Process in groups this large')
    parser.add_argument('--npixel', type=int, default=512, help='Number of pixels in image')
    parser.add_argument('--seed', type=int, default=18051955, help='Random number seed')
    parser.add_argument('--snapshot', type=str, default='False', help='Do snapshot only?')
    parser.add_argument('--opposite', type=str, default='False', help='Move source to opposite side of pointing centre')
    parser.add_argument('--pbradius', type=float, default=4.0, help='Radius of sources to include (in HWHM)')
    parser.add_argument('--pbtype', type=str, default='MID', help='Primary beam model: MID or MID_GAUSS')
    parser.add_argument('--use_agg', type=str, default="True", help='Use Agg matplotlib backend?')
    parser.add_argument('--tsys', type=float, default=0.0, help='System temperature: standard 20K')
    parser.add_argument('--scale', type=float, nargs=2, default=[1.0, 1.0], help='Scale errors by this amount')
    parser.add_argument('--use_radec', type=str, default="False", help='Calculate in RADEC (false)?')
    parser.add_argument('--use_natural', type=str, default="False", help='Use natural weighting?')
    parser.add_argument('--integration_time', type=float, default=600.0, help="Integration time (s)")
    parser.add_argument('--time_range', type=float, nargs=2, default=[-6.0, 6.0], help="Hourangle range (hours")
    parser.add_argument('--time_series', type=str, default='', help="'wind' or 'tracking' or ''")
    parser.add_argument('--pointing_file', type=str, default=None, help="Pointing file")
    parser.add_argument('--time_chunk', type=float, default=1800.0, help="Time for a chunk (s)")
    parser.add_argument('--reference_pointing', type=str, default="False", help="Use reference pointing")
    parser.add_argument('--pointing_directory', type=str, default='../../pointing_error_models/PSD_data/precision/',
                        help='Location of pointing files')
    
    args = parser.parse_args()
    
    use_agg = args.use_agg == "True"
    if use_agg:
        import matplotlib as mpl
        
        mpl.use('Agg')
    from matplotlib import pyplot as plt
    
    use_radec = args.use_radec == "True"
    use_natural = args.use_natural == "True"
    time_series = args.time_series
    pointing_file = args.pointing_file
    scale = numpy.array(args.scale)
    tsys = args.tsys
    integration_time = args.integration_time
    time_range = args.time_range
    time_chunk = args.time_chunk
    snapshot = args.snapshot == 'True'
    opposite = args.opposite == 'True'
    pbtype = args.pbtype
    pbradius = args.pbradius
    reference_pointing = args.reference_pointing == "True"
    pointing_directory = args.pointing_directory
    rmax = args.rmax
    flux_limit = args.flux_limit
    npixel = args.npixel
    
    seed = args.seed
    print("Random number seed is", seed)
    show = args.show == 'True'
    context = args.context
    nworkers = args.nworkers
    nnodes = args.nnodes
    threads_per_worker = args.nthreads
    memory = args.memory
    ngroup = args.ngroup
    
    basename = os.path.basename(os.getcwd())
    
    print("Using %s Dask workers" % nworkers)
    client = get_dask_Client(threads_per_worker=threads_per_worker,
                             processes=threads_per_worker == 1,
                             memory_limit=memory * 1024 * 1024 * 1024,
                             n_workers=nworkers)
    arlexecute.set_client(client=client)
    # n_workers is only relevant if we are using LocalCluster (i.e. a single node) otherwise
    # we need to read the actual number of workers
    actualnworkers = len(arlexecute.client.scheduler_info()['workers'])
    print("Actual number of workers is %d" % actualnworkers)
    nworkers = actualnworkers
    
    time_started = time.time()
    
    # Set up details of simulated observation
    nfreqwin = 1
    diameter = 15.0
    frequency = [1.4e9]
    channel_bandwidth = [1e7]
    
    # Do each 30 minutes in parallel
    start_times = numpy.arange(time_range[0] * 3600, time_range[1] * 3600, time_chunk)
    end_times = start_times + time_chunk
    print("Start times for chunks:")
    pp.pprint(start_times)
    
    times = [numpy.arange(start_times[itime], end_times[itime], integration_time) for itime in range(len(start_times))]
    print("Observation times:")
    pp.pprint(times)
    s2r = numpy.pi / (12.0 * 3600)
    rtimes = s2r * numpy.array(times)
    ntimes = len(rtimes.flat)
    nchunks = len(start_times)
    
    print('%d integrations of duration %.1f s processed in %d chunks' % (ntimes, integration_time, nchunks))
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    mid = create_configuration_from_MIDfile('../../shared/ska1mid_local.cfg', rmax=rmax, location=location)
    
    bvis_graph = [arlexecute.execute(create_blockvisibility)(mid, rtimes[itime], frequency=frequency,
                                                             channel_bandwidth=channel_bandwidth, weight=1.0,
                                                             phasecentre=phasecentre,
                                                             polarisation_frame=PolarisationFrame("stokesI"),
                                                             zerow=True)
                  for itime in range(nchunks)]
    future_bvis_list = arlexecute.persist(bvis_graph, sync=True)
    
    bvis_list0 = arlexecute.compute(bvis_graph[0], sync=True)
    memory_use['bvis_list'] = nchunks * bvis_list0.size()
    del bvis_list0

    vis_graph = [arlexecute.execute(convert_blockvisibility_to_visibility)(bv) for bv in future_bvis_list]
    future_vis_list = arlexecute.persist(vis_graph, sync=True)
    
    vis_list0 = arlexecute.compute(vis_graph[0], sync=True)
    memory_use['vis_list'] = nchunks * vis_list0.size()
    del vis_list0
    
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
    
    FOV_deg = 10.0 * HWHM_deg
    print('%s: HWHM beam = %g deg' % (pbtype, HWHM_deg))
    
    advice_list = arlexecute.execute(advise_wide_field)(future_vis_list[0], guard_band_image=1.0,
                                                        delA=0.02)
    advice = arlexecute.compute(advice_list, sync=True)
    pb_npixel = 1024
    d2r = numpy.pi / 180.0
    pb_cellsize = d2r * FOV_deg / pb_npixel
    cellsize = advice['cellsize']
    
    if show:
        future_vis_list = arlexecute.compute(future_vis_list, sync=True)
        plt.clf()
        for ivis in range(nchunks):
            vis = future_vis_list[ivis]
            plt.plot(-vis.u, -vis.v, '.', color='b', markersize=0.2)
            plt.plot(vis.u, vis.v, '.', color='b', markersize=0.2)
        plt.xlabel('U (wavelengths)')
        plt.ylabel('V (wavelengths)')
        plt.title('UV coverage')
        plt.savefig('uvcoverage.png')
        plt.show(block=False)
        future_vis_list = arlexecute.scatter(future_vis_list)
        
        plt.clf()
        future_bvis_list = arlexecute.compute(future_bvis_list, sync=True)
        for ivis in range(nchunks):
            bvis = future_bvis_list[ivis]
            ha = numpy.pi * bvis.time / 43200.0
            dec = phasecentre.dec.rad
            latitude = bvis.configuration.location.lat.rad
            az, el = hadec_to_azel(ha, dec, latitude)
            plt.plot(ha, az, '.', color='r', label='Azimuth (rad)')
            plt.plot(ha, el, '.', color='b', label='Elevation (rad)')
        plt.xlabel('HA (rad)')
        plt.savefig('azel.png')
        plt.show(block=False)
        future_bvis_list = arlexecute.scatter(future_bvis_list)

    # Construct the skycomponents
    if context == 'singlesource':
        print("Constructing single component")
        offset = [HWHM_deg, 0.0]
        if opposite:
            offset = [-1.0 * offset[0], -1.0 * offset[1]]
        print("Offset from pointing centre = %.3f, %.3f deg" % (offset[0], offset[1]))
        
        # The point source is offset to approximately the halfpower point
        offset_direction = SkyCoord(ra=(+15.0 + offset[0]) * u.deg,
                                    dec=(-45.0 + offset[1]) * u.deg,
                                    frame='icrs', equinox='J2000')
        
        original_components = [Skycomponent(flux=[[1.0]], direction=offset_direction, frequency=frequency,
                                            polarisation_frame=PolarisationFrame('stokesI'))]
        print(original_components[0])
    
    else:
        # Make a skymodel from S3
        print("Constructing s3sky components")
        from wrappers.serial.simulation.testing_support import create_test_skycomponents_from_s3
        
        original_components = create_test_skycomponents_from_s3(flux_limit=flux_limit,
                                                                phasecentre=phasecentre,
                                                                polarisation_frame=PolarisationFrame("stokesI"),
                                                                frequency=numpy.array(frequency),
                                                                radius=pbradius * HWHM)
        print("Created %d original components" % len(original_components))
        # Primary beam points to the phasecentre
        offset_direction = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
        
        
    if time_series == '':
        pes = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
    else:
        pes = ['El15Az0',
               'El15Az135',
               'El15Az180',
               'El15Az45',
               'El15Az90',
               'El45Az0',
               'El45Az135',
               'El45Az180',
               'El45Az45',
               'El45Az90',
               'El90Az0',
               'El90Az135',
               'El90Az180',
               'El90Az45',
               'El90Az90']
        
        
    nants = len(mid.names)
    nbaselines = nants * (nants - 1) // 2

    memory_use['model_list'] = 8 * npixel * npixel * len(frequency) * len(original_components) / 1024 / 1024 / 1024
    memory_use['vp_list'] = 16 * npixel * npixel * len(frequency) * nchunks / 1024 / 1024 / 1024
    print("Memory use (GB)")
    pp.pprint(memory_use)
    total_memory_use = numpy.sum([memory_use[key] for key in memory_use.keys()])

    pp.pprint(arlexecute.client.who_has())
    
    print("Summary of processing:")
    print("    There are %d workers" % nworkers)
    print("    There are %d separate visibility time chunks being processed" % len(future_vis_list))
    print("    The integration time within each chunk is %.1f (s)" % integration_time)
    print("    There are a total of %d integrations" % ntimes)
    print("    There are %d baselines" % nbaselines)
    print("    There are %d components" % len(original_components))
    print("    %d pointing scenario(s) will be tested" % len(pes))
    ntotal = ntimes * nbaselines * len(original_components) * len(pes)
    print("    Total processing %g times-baselines-components-scenarios" % ntotal)
    print("    Approximate total memory use for data = %.3f GB" % total_memory_use)

    # Uniform weighting
    future_model_list = [arlexecute.execute(create_image_from_visibility)(future_vis_list[0], npixel=npixel,
                                                                          frequency=frequency,
                                                                          nchan=nfreqwin, cellsize=cellsize,
                                                                          phasecentre=offset_direction,
                                                                          polarisation_frame=PolarisationFrame(
                                                                              "stokesI"))
                         for i, _ in enumerate(original_components)]
    future_model_list = arlexecute.persist(future_model_list)
    
    psf_list = [arlexecute.execute(create_image_from_visibility)(v, npixel=npixel, frequency=frequency,
                                                                 nchan=nfreqwin, cellsize=cellsize,
                                                                 phasecentre=phasecentre,
                                                                 polarisation_frame=PolarisationFrame("stokesI"))
                for v in future_vis_list]
    psf_list = arlexecute.compute(psf_list, sync=True)
    future_psf_list = arlexecute.scatter(psf_list)
    del psf_list
    
    
    if use_natural:
        print("Using natural weighting")
    else:
        print("Using uniform weighting")
        
        vis_list = weight_list_arlexecute_workflow(future_vis_list, future_psf_list)
        vis_list = arlexecute.compute(vis_list, sync=True)
        future_vis_list = arlexecute.scatter(vis_list)
        del vis_list
        
        bvis_list = [arlexecute.execute(convert_visibility_to_blockvisibility)(vis) for vis in future_vis_list]
        bvis_list = arlexecute.compute(bvis_list, sync=True)
        future_bvis_list = arlexecute.scatter(bvis_list)
        del bvis_list
    
    print("Inverting to get PSF")
    psf_list = invert_list_arlexecute_workflow(future_vis_list, future_psf_list, '2d', dopsf=True)
    psf_list = arlexecute.compute(psf_list, sync=True)
    psf, sumwt = sum_invert_results(psf_list)
    print("PSF sumwt ", sumwt)
    export_image_to_fits(psf, 'PSF_arl.fits')
    if show:
        show_image(psf, cm='gray_r', title='%s PSF' % basename, vmin=-0.01, vmax=0.1)
        plt.savefig('PSF_arl.png')
        plt.show(block=False)
    del psf_list
    
    # ### Calculate the voltage pattern without pointing errors
    vp_list = [arlexecute.execute(create_image_from_visibility)(bv, npixel=pb_npixel, frequency=frequency,
                                                                nchan=nfreqwin, cellsize=pb_cellsize,
                                                                phasecentre=phasecentre,
                                                                override_cellsize=False) for bv in future_bvis_list]
    future_vp_list = arlexecute.persist(vp_list)
    del vp_list
    
    # Optionally show the primary beam, with components if the image is in RADEC coords
    if show:
        pb = arlexecute.execute(create_pb)(future_vp_list[0], pbtype, pointingcentre=phasecentre,
                                           use_local=False)
        pb = arlexecute.compute(pb, sync=True)
        print("Primary beam:", pb)
        if pbtype == 'MID_GRASP':
            show_image(pb, title='%s: primary beam' % basename, vmax=0.01, vmin=0.0)
        else:
            show_image(pb, title='%s: primary beam' % basename, components=original_components, vmax=0.01, vmin=0.0)

        plt.savefig('PB_arl.png')
        export_image_to_fits(pb, 'PB_arl.fits')
        plt.show(block=False)

    # Construct the voltage patterns
    print("Constructing voltage pattern")
    vp_list = [arlexecute.execute(create_vp)(vp, pbtype, pointingcentre=phasecentre, use_local=not use_radec)
               for vp in future_vp_list]
    future_vp_list = arlexecute.persist(vp_list)
    del vp_list
    
    # Make a set of seeds, one per bvis, to ensure that we can get the same errors on different passes
    seeds = numpy.round(numpy.random.uniform(1, numpy.power(2, 31), len(future_bvis_list))).astype(('int'))
    print("Seeds per chunk:")
    pp.pprint(seeds)

    # Set up null pointing error and gain tables. The gaintable will have the voltage pattern in it and so
    # the prediction step will be as if the primary beam had been applied. We need one distinct gain table for each
    # component and for each visibility. Process the data in chunks to avoid needing memory for the entire set of
    # images in one pass
    print("Creating visibilities without any errors")
    print("Predicting error-free visibilities in chunks of %d skymodels" % ngroup)
    dirty_list = list()
    chunk_bvis = [future_bvis_list[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
    chunk_vp_list = [future_vp_list[i:i + ngroup] for i in range(0, len(future_vp_list), ngroup)]
    chunk_seeds = [seeds[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
    for ichunk, chunk in enumerate(chunk_bvis):
        print("Processing chunk %d" % ichunk)
        chunk_dirty_list = create_vis_list_with_errors(chunk_bvis[ichunk], original_components,
                                                       model_list=future_model_list,
                                                       vp_list=chunk_vp_list[ichunk],
                                                       use_radec=use_radec,
                                                       seeds=chunk_seeds[ichunk],
                                                       reference_pointing=reference_pointing)
        result = arlexecute.compute(chunk_dirty_list, sync=True)
        for r in result:
            dirty_list.append(r)
    no_error_dirty, sumwt = sum_invert_results(dirty_list)
    print("Dirty image sumwt ", sumwt)
    print(qa_image(no_error_dirty))
    export_image_to_fits(no_error_dirty, 'dirty_arl.fits')
    if show:
        show_image(no_error_dirty, cm='gray_r', title='%s Dirty image' % basename)  # , vmin=-0.01, vmax=0.1)
        plt.savefig('no_error_dirty_arl.png')
        plt.show(block=False)
    
    results = []
    
    filename = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.csv')
    print('Saving results to %s' % filename)
    plotfile = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.jpg')
    
    epoch = time.strftime("%Y-%m-%d %H:%M:%S")
    
    global_pe = numpy.array(args.global_pe)
    static_pe = args.static_pe
    dynamic_pe = args.dynamic_pe
    
    time_started = time.time()

    # Now loop over all pointing errors
    print("")
    print("***** Starting loop over scenarios ******")
    print("")
    for pe in pes:
        
        result = dict()
        result['context'] = context
        result['nb_name'] = sys.argv[0]
        result['plotfile'] = plotfile
        result['hostname'] = socket.gethostname()
        result['epoch'] = epoch
        result['basename'] = basename
        result['nworkers'] = nworkers
        
        result['npixel'] = npixel
        result['pb_npixel'] = pb_npixel
        result['flux_limit'] = flux_limit
        result['pbtype'] = pbtype
        result['global_pe'] = global_pe
        result['static_pe'] = static_pe
        result['dynamic_pe'] = dynamic_pe
        result['snapshot'] = snapshot
        result['opposite'] = opposite
        result['tsys'] = tsys
        result['scale'] = scale
        result['use_radec'] = use_radec
        result['use_natural'] = use_natural
        result['time_series'] = time_series
        result['integration_time'] = integration_time
        result['seed'] = seed
        result['ntotal'] = ntotal
        
        a2r = numpy.pi / (3600.0 * 180.0)
        if time_series == '':
            global_pointing_error = global_pe
            static_pointing_error = static_pe * pe
            pointing_error = dynamic_pe * pe
            result['static_pointing_error'] = static_pointing_error
            result['dynamic_pointing_error'] = pointing_error
            result['global_pointing_error'] = global_pointing_error
        
            print("Pointing errors: global (%.1f, %.1f) arcsec, static %.1f arcsec, dynamic %.1f arcsec" %
                  (global_pointing_error[0], global_pointing_error[1], static_pointing_error,
                   pointing_error))
            file_name = 'PE_%.1f_arcsec_arl' % pe

            error_dirty_list = list()
            chunk_bvis = [future_bvis_list[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
            chunk_vp_list = [future_vp_list[i:i + ngroup] for i in range(0, len(future_vp_list), ngroup)]
            chunk_seeds = [seeds[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
            for ichunk, chunk in enumerate(chunk_bvis):
                print("Processing chunk %d" % ichunk)
                chunk_dirty_list = create_vis_list_with_errors(chunk_bvis[ichunk], original_components,
                                                               model_list=future_model_list,
                                                               vp_list=chunk_vp_list[ichunk],
                                                               use_radec=use_radec,
                                                               pointing_error=a2r * pointing_error,
                                                               static_pointing_error=a2r * static_pointing_error,
                                                               global_pointing_error=a2r * global_pointing_error,
                                                               seeds=chunk_seeds[ichunk])
                this_result = arlexecute.compute(chunk_dirty_list, sync=True)
                for r in this_result:
                    error_dirty_list.append(r)
        else:
            result['pointing_file'] = "%s/%s.dat" % (pointing_directory, pe)
            az = float(pe.split('dat')[0].split('Az')[1])
            el = float(pe.split('dat')[0].split('Az')[0].split('El')[1])
            result['wind_azimuth'] = az
            result['wind_elevation'] = el
            file_name = 'PE_%s_%s_arl' % (pe, time_series)
            print("Pointing errors: type of time series %s from file %s" % (time_series, result['pointing_file']))

            pp.pprint(arlexecute.client.who_has())

            error_dirty_list = list()
            chunk_bvis = [future_bvis_list[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
            chunk_vp_list = [future_vp_list[i:i + ngroup] for i in range(0, len(future_vp_list), ngroup)]
            chunk_seeds = [seeds[i:i + ngroup] for i in range(0, len(future_bvis_list), ngroup)]
            for ichunk, chunk in enumerate(chunk_bvis):
                print("Processing chunk %d" % ichunk)
                chunk_dirty_list = create_vis_list_with_errors(chunk_bvis[ichunk], original_components,
                                                               model_list=future_model_list,
                                                               vp_list=chunk_vp_list[ichunk],
                                                               use_radec=use_radec,
                                                               time_series=time_series,
                                                               pointing_file=result['pointing_file'],
                                                               seeds=chunk_seeds[ichunk],
                                                               reference_pointing=reference_pointing)
                this_result = arlexecute.compute(chunk_dirty_list, sync=True)
                for r in this_result:
                    error_dirty_list.append(r)

        error_dirty, sumwt = sum_invert_results(error_dirty_list)
        print("Dirty image sumwt", sumwt)
        del error_dirty_list
        error_dirty.data -= no_error_dirty.data
        print(qa_image(error_dirty))
        
        if show:
            show_image(error_dirty, cm='gray_r', title='Pointing error %s' % file_name)
            plt.savefig('%s.png' % file_name)
            plt.show(block=False)
        
        qa = qa_image(error_dirty)
        _, _, ny, nx = error_dirty.shape
        for field in ['maxabs', 'rms', 'medianabs']:
            result["onsource_" + field] = qa.data[field]
        result['onsource_abscentral'] = numpy.abs(error_dirty.data[0, 0, ny // 2, nx // 2])
        
        result['elapsed_time'] = time.time() - time_started
        print('Elapsed time = %.1f (s)' % result['elapsed_time'])
        
        results.append(result)
    
    pp.pprint(results)
    
    print("Total processing %g times-baselines-components-scenarios" % ntotal)
    processing_rate = ntotal / (nworkers * (time.time() - time_started))
    # Typical values:
    # Tim-MBP, MacBookPro14,3 Intel Core i7 2.9 GHz, 5818.72 /s/worker
    # Sheldon, Intel(R) Core(TM) i7-6900K CPU @ 3.20GHz, 22000.0 /s/worker
    # CSD3, single node, Intel(R) Xeon(R) Gold 6142 CPU @ 2.60GHz, 29522.8 /s/worker
    # CSD3, multinode, Intel(R) Xeon(R) Gold 6142 CPU @ 2.60GHz, 12600.0 /s/worker
    #
    print("Processing rate of time-baseline-component-scenario = %g per worker-second" % processing_rate)

    with open(filename, 'a') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=results[0].keys(), delimiter=',', quotechar='|',
                                quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for result in results:
            result["processing_rate"] = processing_rate
            writer.writerow(result)
        csvfile.close()
    
    if time_series == '':
        title = '%s, %.3f GHz, %d times: dynamic %g, static %g \n%s %s %s' % \
                (context, frequency[0] * 1e-9, ntimes, dynamic_pe, static_pe, socket.gethostname(), epoch,
                 basename)
        plt.clf()
        colors = ['b', 'r', 'g', 'y']
        for ifield, field in enumerate(['onsource_maxabs', 'onsource_rms', 'onsource_medianabs']):
            plt.loglog(pes, [1e6 * result[field] for result in results], '-', label=field, color=colors[ifield])
        
        plt.xlabel('Pointing error (arcsec)')
        plt.ylabel('Error (uJy)')
        
        plt.title(title)
        plt.legend(fontsize='x-small')
        print('Saving plot to %s' % plotfile)
        
        plt.savefig(plotfile)
        plt.show(block=False)
        
    else:
    
        title = '%s, %.3f GHz, %d times %s \n%s %s %s' % \
                (context, frequency[0] * 1e-9, ntimes, time_series, socket.gethostname(), epoch,
                 basename)
        bar_width = 0.35
        opacity = 0.8

        plt.clf()
        index = numpy.arange(len(pes))
        fig, ax = plt.subplots()
        colors = ['b', 'r', 'g', 'y']
        for ifield, field in enumerate(['onsource_rms', 'onsource_medianabs']):

            plt.bar(index + ifield * bar_width, [1e6 * result[field] for result in results],
                    bar_width, label=field, color=colors[ifield],
                    alpha=opacity)
    
        plt.xlabel('Pointing file')
        plt.ylabel('Error (uJy)')
        plt.xticks(numpy.arange(len(pes))+ bar_width, pes, rotation='vertical')
        plt.title(title)
        plt.legend(fontsize='x-small')
        print('Saving plot to %s' % plotfile)
    
        plt.tight_layout()
        plt.savefig(plotfile)
        plt.show(block=False)

