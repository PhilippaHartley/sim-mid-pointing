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

# mpl.use('Agg')

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
    
    args = parser.parse_args()
    
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
    ntimes = 65
    diameter = 15.0
    frequency = [1.4e9]
    channel_bandwidth = [1e7]
    
    h2r = numpy.pi / 12.0
    times = numpy.linspace(-6 * h2r, +6 * h2r, ntimes)
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    outlier_phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-35.0 * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    mid = create_configuration_from_MIDfile('../../shared/ska1mid.cfg', rmax=rmax, location=location)

    # ARL calculates uvw incorrectly so we use those from CASA
    from processing_components.visibility.base import create_blockvisibility_from_ms
    meqbvis = create_blockvisibility_from_ms('../../shared/PE_0.0_arcsec.ms')[0]
    meqbvis.phasecentre = phasecentre
    meqvis = convert_blockvisibility_to_visibility(meqbvis)
    print(meqbvis)

    block_vis = create_blockvisibility(mid, times, frequency=frequency,
                                           channel_bandwidth=channel_bandwidth, weight=1.0,
                                           phasecentre=phasecentre,
                                           polarisation_frame=PolarisationFrame("stokesI"), zerow=False)
    print(block_vis)
    vis = convert_blockvisibility_to_visibility(block_vis)
    
    print(numpy.max(numpy.abs(vis.data['uvw'] - meqvis.uvw)))
    vis.data['uvw'] = meqvis.uvw
    print(numpy.max(numpy.abs(vis.data['uvw'] - meqvis.uvw)))

    print(numpy.max(numpy.abs(vis.data['antenna1']-meqvis.antenna1)))
    print(numpy.max(numpy.abs(vis.data['antenna2']-meqvis.antenna2)))

    vis.data['antenna1'] = meqvis.antenna1
    vis.data['antenna2'] = meqvis.antenna2
    vis.data['uvw'][...,2] = 0.0
    
    advice = advise_wide_field(vis, guard_band_image=1.0, delA=0.02)
    
    # We need the HWHM of the primary beam. Got this by trial and error
    HWHM_deg = 1.03 * 180.0 * 3e8 / (numpy.pi * diameter * frequency[0])
    
    print('HWHM beam = %g deg' % HWHM_deg)
    HWHM = HWHM_deg * numpy.pi / 180.0
    
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
        show_image(psf, cm='gray_r', title='PSF')
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
        plt.show(block=False)
    
    print("Constructing voltage pattern")
    vp = create_vp(vp, 'MID', pointingcentre=pb_direction)
    pt = create_pointingtable_from_blockvisibility(block_vis, vp)
    
    no_error_pt = simulate_pointingtable(pt, 0.0, 0.0, seed=18051955)
    export_pointingtable_to_hdf5(no_error_pt, 'pointingsim_%s_noerror_pointingtable.hdf5' % context)
    no_error_gt = create_gaintable_from_pointingtable(block_vis, original_components, no_error_pt, vp)
    
    no_error_sm = [SkyModel(components=[original_components[i]], gaintable=no_error_gt[i])
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
    
    pes = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
    results = []
    
    filename = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.csv')
    print('Saving results to %s' % filename)
    plotfile = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.jpg')
    
    epoch = time.strftime("%Y-%m-%d %H:%M:%S")
    
    global_pe = numpy.array(args.global_pe)
    static_pe = args.static_pe
    dynamic_pe = args.dynamic_pe
    
    # Now loop over all pointing errors
    for pe in pes:
        
        result = dict()
        result['context'] = context
        result['nb_name'] = sys.argv[0]
        result['plotfile'] = plotfile
        result['hostname'] = socket.gethostname()
        result['epoch'] = epoch
        result['npixel'] = npixel
        result['pb_npixel'] = pb_npixel
        result['flux_limit'] = flux_limit
        
        result['global_pe'] = global_pe
        result['static_pe'] = static_pe
        result['dynamic_pe'] = dynamic_pe
        
        a2r = numpy.pi / (3600.0 * 180.0)
        global_pointing_error = global_pe
        static_pointing_error = static_pe * pe
        pointing_error = dynamic_pe * pe
        
        result['static_pointing_error'] = static_pointing_error
        result['dynamic_pointing_error'] = pointing_error
        result['global_pointing_error'] = global_pointing_error
        
        print("Pointing errors: global (%.1f, %.1f) arcsec, static %.1f arcsec, dynamic %.1f arcsec" %
              (global_pointing_error[0], global_pointing_error[1], static_pointing_error,
               pointing_error))
        
        error_pt = simulate_pointingtable(pt, pointing_error=pointing_error * a2r,
                                          static_pointing_error=static_pointing_error * a2r,
                                          global_pointing_error=global_pointing_error * a2r, seed=18051955)
        export_pointingtable_to_hdf5(error_pt,
                                     'pointingsim_%s_error_%.0farcsec_pointingtable.hdf5' % (context, pe))
        
        error_gt = create_gaintable_from_pointingtable(block_vis, original_components, error_pt, vp)
        
        error_sm = [SkyModel(components=[original_components[i]], gaintable=error_gt[i])
                    for i, _ in enumerate(original_components)]
        
        error_blockvis = copy_visibility(block_vis, zero=True)
        
        print("Predicting corrupted visibilities in chunks of %d skymodels" % ngroup)
        future_vis = arlexecute.scatter(error_blockvis)
        chunks = [error_sm[i:i + ngroup] for i in range(0, len(error_sm), ngroup)]
        for chunk in chunks:
            temp_vis = predict_skymodel_list_compsonly_arlexecute_workflow(future_vis, chunk, context='2d', docal=True)
            work_vis = arlexecute.compute(temp_vis, sync=True)
            for w in work_vis:
                error_blockvis.data['vis'] += w.data['vis']
            assert numpy.max(numpy.abs(error_blockvis.data['vis'])) > 0.0
        
        error_blockvis.data['vis'] -= no_error_blockvis.data['vis']
        
        print("Inverting to get on-source dirty image")
        error_vis = convert_blockvisibility_to_visibility(error_blockvis)
        dirty = invert_list_arlexecute_workflow([error_vis], [model], '2d')
        dirty, sumwt = arlexecute.compute(dirty, sync=True)[0]
        
        export_image_to_fits(dirty, 'PE_%.1f_arcsec_arl.fits' % pe)
        
        if show:
            show_image(dirty, cm='gray_r', title='Pointing error %.1f arcsec ARL' % pe)
            plt.savefig('PE_%.1f_arcsec_arl.png' % pe)
            plt.show(block=False)
        
        qa = qa_image(dirty)
        _, _, ny, nx = dirty.shape
        for field in ['maxabs', 'rms', 'medianabs']:
            result["onsource_" + field] = qa.data[field]
        result['onsource_abscentral'] = numpy.abs(dirty.data[0, 0, ny // 2, nx // 2])
        
        if outlier:
            print("Inverting to get outlier dirty image")
            outlier_model = create_image_from_visibility(error_vis, npixel=npixel, frequency=frequency,
                                                         nchan=nfreqwin, cellsize=cellsize,
                                                         phasecentre=outlier_phasecentre)
            outlier_dirty = invert_list_arlexecute_workflow([error_vis], [outlier_model], '2d')
            outlier_dirty, outlier_sumwt = arlexecute.compute(outlier_dirty, sync=True)[0]
            
            if show:
                show_image(outlier_dirty, cm='gray_r', title='Outlier residual image (dec -35deg)')
                plt.show(block=False)
            
            qa = qa_image(outlier_dirty)
            _, _, ny, nx = outlier_dirty.shape
            for field in ['maxabs', 'rms', 'medianabs']:
                result["offsource_" + field] = qa.data[field]
            result['offsource_abscentral'] = numpy.abs(outlier_dirty.data[0, 0, ny // 2, nx // 2])
        
        results.append(result)
    
    import pprint
    
    pp = pprint.PrettyPrinter()
    pp.pprint(results)
    
    with open(filename, 'a') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=results[0].keys(), delimiter=',', quotechar='|',
                                quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for result in results:
            writer.writerow(result)
        csvfile.close()
    
    plt.clf()
    colors = ['b', 'r', 'g', 'y']
    for ifield, field in enumerate(['onsource_maxabs', 'onsource_rms', 'onsource_medianabs', 'onsource_abscentral']):
        plt.loglog(pes, [result[field] for result in results], '-', label=field, color=colors[ifield])
    
    if outlier:
        for ifield, field in enumerate(['outlier_maxabs', 'outlier_rms', 'outlier_medianabs', 'offsource_abscentral']):
            plt.loglog(pes, [result[field] for result in results], '--', label=field, color=colors[ifield])
        plt.xlabel('Pointing error (arcsec)')
        plt.ylabel('Error (Jy)')
    
    title = '%s, %.3f GHz, %d times: dynamic %g, static %g \n%s %s' % \
            (context, frequency[0] * 1e-9, ntimes, dynamic_pe, static_pe, socket.gethostname(), epoch)
    plt.title(title)
    plt.legend(fontsize='x-small')
    print('Saving plot to %s' % plotfile)
    
    plt.savefig(plotfile)
    plt.show()
    plt.show()
