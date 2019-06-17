"""Simulation of the effect of pointing errors on MID observations
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

from wrappers.serial.visibility.base import create_blockvisibility
from wrappers.serial.image.operations import show_image, qa_image, export_image_to_fits
from wrappers.serial.simulation.configurations import create_configuration_from_MIDfile
from wrappers.serial.simulation.testing_support import simulate_pointingtable, simulate_pointingtable_from_timeseries
from wrappers.serial.simulation.noise import addnoise_visibility
from wrappers.serial.imaging.primary_beams import create_vp, create_pb
from wrappers.serial.imaging.base import create_image_from_visibility, advise_wide_field
from wrappers.serial.calibration.pointing import create_pointingtable_from_blockvisibility
from wrappers.serial.simulation.pointing import simulate_gaintable_from_pointingtable
from wrappers.arlexecute.visibility.base import copy_visibility
from wrappers.arlexecute.visibility.coalesce import convert_blockvisibility_to_visibility, \
    convert_visibility_to_blockvisibility
from processing_library.util.coordinate_support import hadec_to_azel

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

def create_vis_list_with_errors(bvis_list, original_components, use_radec=False, pointing_error=0.0,
                                static_pointing_error=0.0,
                                global_pointing_error=None, time_series=''):
    if global_pointing_error is None:
        global_pointing_error = [0.0, 0.0]
    
    # One pointing table per visibility
    pt_list = [arlexecute.execute(create_pointingtable_from_blockvisibility)(bvis) for bvis in bvis_list]
    
    if time_series is '':
        error_pt_list = [arlexecute.execute(simulate_pointingtable)(pt, pointing_error=pointing_error,
                                                                    static_pointing_error=static_pointing_error,
                                                                    global_pointing_error=global_pointing_error)
                         for pt in pt_list]
    else:
        error_pt_list = [arlexecute.execute(simulate_pointingtable_from_timeseries)(pt, type=time_series)
                         for pt in pt_list]
    
    if show:
        error_pt_list = arlexecute.compute(error_pt_list, sync=True)
        plt.clf()
        r2a = 180.0 * 3600.0 / numpy.pi
        for pt in error_pt_list:
            plt.plot(pt.time, r2a * pt.pointing[:, 0, 0, 0, 0], '.', color='r')
            plt.plot(pt.time, r2a * pt.pointing[:, 0, 0, 0, 1], '.', color='b')
        plt.title("%s: pointing" % basename)
        plt.xlabel('Time (s)')
        plt.ylabel('Offset (arcsec)')
        plt.show()
    
    # Create the gain tables, one per Visibility and per component
    error_gt_list = [arlexecute.execute(simulate_gaintable_from_pointingtable)
                     (bvis, original_components, error_pt_list[ibv], vp_list[ibv], use_radec=use_radec)
                     for ibv, bvis in enumerate(bvis_list)]
    if show:
        error_gt_list = arlexecute.compute(error_gt_list, sync=True)
        plt.clf()
        for gt in error_gt_list:
            plt.plot(gt[0].time, 1.0 / numpy.abs(gt[0].gain[:, 0, 0, 0]), '.')
        plt.ylim([0.0, 1.1])
        plt.title("%s: amplitude gain" % basename)
        plt.xlabel('Time (s)')
        plt.show()
    
    # Each component in original components becomes a separate skymodel
    # Inner nest is over skymodels, outer is over bvis's
    error_sm_list = [[
        arlexecute.execute(SkyModel, nout=1)(components=[original_components[i]], gaintable=error_gt_list[ibv][i])
        for i, _ in enumerate(original_components)] for ibv, bv in enumerate(bvis_list)]
    
    # Predict each visibility for each skymodel. We keep all the visibilities separate
    # and add up dirty images at the end of processing. We calibrate which applies the voltage pattern
    error_bvis_list = [arlexecute.execute(copy_visibility, nout=1)(bvis, zero=True) for bvis in bvis_list]
    error_bvis_list = [predict_skymodel_list_compsonly_arlexecute_workflow(error_bvis_list[ibv], error_sm_list[ibv],
                                                                           context='2d', docal=True)
                       for ibv, bvis in enumerate(error_bvis_list)]
    # Inner nest is over skymodels, outer is over bvis's
    return error_bvis_list


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
    parser.add_argument('--seed', type=int, default=18051955, help='Random number seed')
    parser.add_argument('--snapshot', type=str, default='False', help='Do snapshot only?')
    parser.add_argument('--opposite', type=str, default='False', help='Move source to opposite side of pointing centre')
    parser.add_argument('--pbtype', type=str, default='MID', help='Primary beam model: MID or MID_GAUSS')
    parser.add_argument('--use_agg', type=str, default="True", help='Use Agg matplotlib backend?')
    parser.add_argument('--tsys', type=float, default=0.0, help='System temperature: standard 20K')
    parser.add_argument('--scale', type=float, nargs=2, default=[1.0, 1.0], help='Scale errors by this amount')
    parser.add_argument('--use_radec', type=str, default="False", help='Calculate in RADEC (false)?')
    parser.add_argument('--use_natural', type=str, default="False", help='Use natural weighting?')
    parser.add_argument('--integration_time', type=float, default=600.0, help="Integration time (s)")
    parser.add_argument('--time_range', type=float, nargs=2, default=[-6.0, 6.0], help="Hourangle range (hours")
    parser.add_argument('--time_series', type=str, default='', help="'wind' or 'tracking' or ''")
    
    args = parser.parse_args()
    
    scale = numpy.array(args.scale)
    tsys = args.tsys
    integration_time = args.integration_time
    time_range = args.time_range
    
    use_agg = args.use_agg == "True"
    use_radec = args.use_radec == "True"
    use_natural = args.use_natural == "True"
    time_series = args.time_series
    
    if use_agg:
        import matplotlib as mpl
        
        mpl.use('Agg')
    
    from matplotlib import pyplot as plt
    
    snapshot = args.snapshot == 'True'
    opposite = args.opposite == 'True'
    pbtype = args.pbtype
    
    seed = args.seed
    
    print("Random number seed is", seed)
    numpy.random.seed(seed)
    show = args.show == 'True'
    context = args.context
    rmax = args.rmax
    flux_limit = args.flux_limit
    
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
    
    # Set up details of simulated observation
    nfreqwin = 1
    diameter = 15.0
    frequency = [1.4e9]
    channel_bandwidth = [1e7]
    
    # Do each 30 minutes in parallel
    start_times = numpy.arange(time_range[0] * 3600, time_range[1] * 3600, 1800.0)
    end_times = start_times + 1800.0
    print(start_times)
    
    times = [numpy.arange(start_times[itime], end_times[itime], integration_time) for itime in range(len(start_times))]
    print(times)
    s2r = numpy.pi / (12.0 * 3600)
    rtimes = s2r * numpy.array(times)
    ntimes = len(rtimes.flat)
    nchunks = len(start_times)

    print('%d integrations processed in %d chunks' % (ntimes, nchunks))
    
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    mid = create_configuration_from_MIDfile('../../shared/ska1mid_local.cfg', rmax=rmax, location=location)
    
    bvis_list = [arlexecute.execute(create_blockvisibility)(mid, rtimes[itime], frequency=frequency,
                                                            channel_bandwidth=channel_bandwidth, weight=1.0,
                                                            phasecentre=phasecentre,
                                                            polarisation_frame=PolarisationFrame("stokesI"), zerow=True)
                 for itime in range(nchunks)]
    bvis_list = arlexecute.compute(bvis_list, sync=True)
    vis_list = [arlexecute.execute(convert_blockvisibility_to_visibility)(bv) for bv in bvis_list]
    vis_list = arlexecute.persist(vis_list, sync=True)
    
    # We need the HWHM of the primary beam. Got this by trial and error
    if pbtype == 'MID':
        HWHM_deg = 0.6 * 1.4e9 / frequency[0]
    elif pbtype == 'MID_GRASP':
        HWHM_deg = 0.755 * 1.4e9 / frequency[0]
    elif pbtype == 'MID_GAUSS':
        HWHM_deg = 0.6 * 1.4e9 / frequency[0]
    else:
        HWHM_deg = 0.6 * 1.4e9 / frequency[0]
    
    FOV_deg = 10.0 * HWHM_deg
    print('%s: HWHM beam = %g deg' % (pbtype, HWHM_deg))
    
    advice_list = [arlexecute.execute(advise_wide_field)(v, guard_band_image=1.0, delA=0.02) for v in vis_list]
    advice = arlexecute.compute(advice_list, sync=True)[0]
    pb_npixel = 1024
    d2r = numpy.pi / 180.0
    pb_cellsize = d2r * FOV_deg / pb_npixel
    cellsize = advice['cellsize']
    npixel = 512
    
    pp.pprint(client.has_what())
    if False:
        plt.clf()
        for vis in vis_list:
            plt.plot(-vis.u, -vis.v, '.', color='b', markersize=0.2)
            plt.plot(vis.u, vis.v, '.', color='b', markersize=0.2)
        plt.xlabel('U (wavelengths)')
        plt.ylabel('V (wavelengths)')
        plt.title('UV coverage')
        plt.savefig('uvcoverage.png')
        plt.show(block=False)
        
        plt.clf()
        for bvis in bvis_list:
            ha = numpy.pi * bvis.time / 43200.0
            dec = phasecentre.dec.rad
            latitude = bvis.configuration.location.lat.rad
            az, el = hadec_to_azel(ha, dec, latitude)
            plt.plot(ha, az, '.', color='r', label='Azimuth (rad)')
            plt.plot(ha, el, '.', color='b', label='Elevation (rad)')
        plt.xlabel('HA (rad)')
        plt.savefig('azel.png')
        plt.show(block=False)
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
                                                                radius=pb_cellsize * pb_npixel / 2.0)
        # Primary beam points to the phasecentre
        offset_direction = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    
    # Uniform weighting
    pp.pprint(client.has_what())


    psf_list = [arlexecute.execute(create_image_from_visibility)(v, npixel=npixel, frequency=frequency,
                                                                 nchan=nfreqwin, cellsize=cellsize,
                                                                 phasecentre=phasecentre,
                                                                 polarisation_frame=PolarisationFrame("stokesI"))
                for v in vis_list]
    psf_list = arlexecute.compute(psf_list, sync=True)
    
    if use_natural:
        print("Using natural weighting")
    else:
        print("Using uniform weighting")
        pp.pprint(client.has_what())

        vis_list = weight_list_arlexecute_workflow(vis_list, psf_list)
        vis_list = arlexecute.compute(vis_list, sync=True)
        bvis_list = [arlexecute.execute(convert_visibility_to_blockvisibility)(vis) for vis in vis_list]
        bvis_list = arlexecute.compute(bvis_list, sync=True)
    
    print("Inverting to get on-source PSF")
    psf_list = invert_list_arlexecute_workflow(vis_list, psf_list, '2d', dopsf=True)
    psf_list = arlexecute.compute(psf_list, sync=True)
    psf, sumwt = sum_invert_results(psf_list)
    print("PSF sumwt ", sumwt)
    export_image_to_fits(psf, 'PSF_arl.fits')
    if show:
        show_image(psf, cm='gray_r', title='%s PSF' % basename, vmin=-0.01, vmax=0.1)
        plt.savefig('PSF_arl.png')
        plt.show(block=False)
    
    model_list = [arlexecute.execute(create_image_from_visibility)(v, npixel=npixel, frequency=frequency,
                                                                   nchan=nfreqwin, cellsize=cellsize,
                                                                   phasecentre=offset_direction,
                                                                   polarisation_frame=PolarisationFrame("stokesI"))
                  for v in vis_list]
    
    # ### Calculate the voltage pattern without pointing errors
    vp_list = [arlexecute.execute(create_image_from_visibility)(bv, npixel=pb_npixel, frequency=frequency,
                                                                nchan=nfreqwin, cellsize=pb_cellsize,
                                                                phasecentre=phasecentre,
                                                                override_cellsize=False) for bv in bvis_list]
    vp_list = arlexecute.compute(vp_list, sync=True)
    
    if show:
        pb = create_pb(vp_list[0], pbtype, pointingcentre=phasecentre, use_local=not use_radec)
        print("Primary beam:", pb)
        show_image(pb, title='%s: primary beam' % basename)
        plt.savefig('PB_arl.png')
        export_image_to_fits(pb, 'PB_arl.fits')
        plt.show(block=False)
    
    print("Constructing voltage pattern")
    vp_list = [arlexecute.execute(create_vp)(vp, pbtype, pointingcentre=phasecentre, use_local=not use_radec)
               for vp in vp_list]
    vp_list = arlexecute.compute(vp_list, sync=True)
    print("Voltage pattern:", vp_list[0])
    
    # Set up null pointing error and gain tables. The gaintable will have the voltage pattern in it and so
    # the prediction step will be as if the primary beam had been applied. We need one distinct gain table for each
    # component and for each visibility.
    
    no_error_bvis_list = create_vis_list_with_errors(bvis_list, original_components, use_radec=use_radec)
    no_error_bvis_list = arlexecute.compute(no_error_bvis_list, sync=True)

    no_error_vis_list = [[arlexecute.execute(convert_blockvisibility_to_visibility)(no_error_bvis_list[inebvl][ineb])
                       for ineb, _ in enumerate(nebvl)] for inebvl, nebvl in enumerate(no_error_bvis_list)]
    no_error_vis_list = arlexecute.compute(no_error_vis_list, sync=True)

    print("After creating no_error_vis_list:")
    pp.pprint(client.has_what())

    # Inner nest is over skymodels, outer is over bvis's: need to swap these
    skymodel0_vis_list = [nevl[0] for nevl in no_error_vis_list]
    print("Inverting to get dirty image")
    dirty_list = invert_list_arlexecute_workflow(skymodel0_vis_list, model_list, '2d')
    dirty_list = arlexecute.compute(dirty_list, sync=True)
    dirty, sumwt = sum_invert_results(dirty_list)
    print("Dirty image sumwt ", sumwt)
    print(qa_image(dirty))
    export_image_to_fits(dirty, 'dirty_arl.fits')
    if show:
        show_image(dirty, cm='gray_r', title='%s Dirty image' % basename)  # , vmin=-0.01, vmax=0.1)
        plt.savefig('dirty_arl.png')
        plt.show(block=False)

    if time_series == '':
        pes = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
    else:
        pes = [0.0]
    
    results = []
    
    filename = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.csv')
    print('Saving results to %s' % filename)
    plotfile = seqfile.findNextFile(prefix='pointingsimulation_%s_' % socket.gethostname(), suffix='.jpg')
    
    epoch = time.strftime("%Y-%m-%d %H:%M:%S")
    
    global_pe = numpy.array(args.global_pe)
    static_pe = args.static_pe
    dynamic_pe = args.dynamic_pe
    
    # Now loop over all pointing errors
    print("***** Starting loop over pointing error ******")
    for pe in pes:
        # We want all simulations to have the same errors, just scaled
        numpy.random.seed(seed)
        
        result = dict()
        result['context'] = context
        result['nb_name'] = sys.argv[0]
        result['plotfile'] = plotfile
        result['hostname'] = socket.gethostname()
        result['epoch'] = epoch
        result['basename'] = basename
        
        result['npixel'] = npixel
        result['pb_npixel'] = pb_npixel
        result['flux_limit'] = flux_limit
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
        result['seed'] = seed
        
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
        
        error_bvis_list = create_vis_list_with_errors(bvis_list, original_components, use_radec=use_radec,
                                                      pointing_error=a2r * pointing_error,
                                                      static_pointing_error=a2r * static_pointing_error,
                                                      global_pointing_error=a2r * global_pointing_error,
                                                      time_series=time_series)
        error_bvis_list = arlexecute.compute(error_bvis_list, sync=True)

        print("After creating error_vis_list:")
        pp.pprint(client.has_what())

        # Dask function to be executed: difference, add noise, and convert to visibility
        def difference_add_noise_convert(error_bvis, no_error_bvis):
            error_bvis.data['vis'] -= no_error_bvis.data['vis']
            error_bvis = addnoise_visibility(error_bvis, tsys)
            error_vis = convert_blockvisibility_to_visibility(error_bvis)
            return error_vis
        
        error_vis_list = [[arlexecute.execute(difference_add_noise_convert)
                           (error_bvis_list[iebvl][ieb], no_error_bvis_list[iebvl][ieb])
                           for ieb, _ in enumerate(ebvl)] for iebvl, ebvl in enumerate(error_bvis_list)]

        error_vis_list = arlexecute.compute(error_vis_list, sync=True)
        
        print("Inverting to get on-source dirty image")
        skymodel0_vis_list = [evl[0] for evl in error_vis_list]
        dirty_list = invert_list_arlexecute_workflow(skymodel0_vis_list, model_list, '2d')
        dirty_list = arlexecute.compute(dirty_list, sync=True)
        dirty, sumwt = sum_invert_results(dirty_list)
        print(qa_image(dirty))
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
        
        results.append(result)
    
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
    
    plt.xlabel('Pointing error (arcsec)')
    plt.ylabel('Error (Jy)')
    
    title = '%s, %.3f GHz, %d times: dynamic %g, static %g \n%s %s %s' % \
            (context, frequency[0] * 1e-9, ntimes, dynamic_pe, static_pe, socket.gethostname(), epoch,
             basename)
    plt.title(title)
    plt.legend(fontsize='x-small')
    print('Saving plot to %s' % plotfile)
    
    plt.savefig(plotfile)
    plt.show()
