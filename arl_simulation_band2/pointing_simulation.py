"""Simulation of the effect of pointing errors on MID observations
"""
import csv
import os
import socket
import sys
import time

import seqfile

from rascil.data_models.parameters import rascil_path

results_dir = rascil_path('test_results')

import numpy

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from rascil.data_models.polarisation import PolarisationFrame
from rascil.data_models.memory_data_models import Skycomponent, SkyModel
from rascil.data_models.data_model_helpers import export_pointingtable_to_hdf5

from rascil.processing_components.visibility.base import create_blockvisibility
from rascil.processing_components.image.operations import show_image, qa_image, export_image_to_fits
from rascil.processing_components.simulation.configurations import create_configuration_from_MIDfile
from rascil.processing_components.simulation.testing_support import simulate_pointingtable, simulate_pointingtable_from_timeseries
from rascil.processing_components.simulation.noise import addnoise_visibility
from rascil.processing_components.imaging.primary_beams import create_vp, create_pb
from rascil.processing_components.imaging.base import create_image_from_visibility, advise_wide_field
from rascil.processing_components.calibration.pointing import create_pointingtable_from_blockvisibility
from rascil.processing_components.simulation.pointing import simulate_gaintable_from_pointingtable
from workflows.rsexecute.visibility.base import copy_visibility
from workflows.rsexecute.visibility.coalesce import convert_blockvisibility_to_visibility, \
    convert_visibility_to_blockvisibility
from processing_library.util.coordinate_support import hadec_to_azel

from rascil.workflows.rsexecute.skymodel.skymodel_rsexecute import predict_skymodel_list_compsonly_rsexecute_workflow
from rascil.workflows.rsexecute.imaging.imaging_rsexecute import invert_list_rsexecute_workflow
from rascil.workflows.serial.imaging.imaging_serial import weight_list_serial_workflow

from workflows.rsexecute.execution_support.rsexecute import rsexecute
from workflows.rsexecute.execution_support.dask_init import get_dask_client

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
    parser.add_argument('--snapshot', type=str, default='False', help='Do snapshot only?')
    parser.add_argument('--opposite', type=str, default='False', help='Move source to opposite side of pointing centre')
    parser.add_argument('--pbtype', type=str, default='MID', help='Primary beam model: MID or MID_GAUSS')
    parser.add_argument('--use_agg', type=str, default="True", help='Use Agg matplotlib backend?')
    parser.add_argument('--tsys', type=float, default=0.0, help='System temperature: standard 20K')
    parser.add_argument('--scale', type=float, nargs=2, default=[1.0, 1.0], help='Scale errors by this amount')
    parser.add_argument('--use_radec', type=str, default="False", help='Calculate in RADEC (false)?')
    parser.add_argument('--use_natural', type=str, default="False", help='Use natural weighting?')
    parser.add_argument('--integration_time', type=float, default=43200.0/65.0, help="Integration time (s)")
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
    
    print("Random number seed is ", seed)
    numpy.random.seed(seed)
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
    
    basename = os.path.basename(os.getcwd())
    
    print("Using %s Dask workers" % nworkers)
    client = get_dask_client(threads_per_worker=threads_per_worker,
                             processes=threads_per_worker == 1,
                             memory_limit=memory * 1024 * 1024 * 1024,
                             n_workers=nworkers)
    rsexecute.set_client(client=client)
    
    # Set up details of simulated observation
    nfreqwin = 1
    ntimes = 65
    diameter = 15.0
    frequency = [1.4e9]
    channel_bandwidth = [1e7]
    
    if snapshot:
        ntimes = 1
        times = [0.0]
    else:
        ntimes = numpy.round(3600.0 * (time_range[1] - time_range[0]) /integration_time).astype('int')
        h2r = numpy.pi / 12.0
        times = numpy.linspace(time_range[0], time_range[1], ntimes)
        
    print("Observations at ", times)
    times *= h2r
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    outlier_phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-35.0 * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    mid = create_configuration_from_MIDfile('../../shared/ska1mid_local.cfg', rmax=rmax, location=location)
    
    block_vis = create_blockvisibility(mid, times, frequency=frequency,
                                       channel_bandwidth=channel_bandwidth, weight=1.0,
                                       phasecentre=phasecentre,
                                       polarisation_frame=PolarisationFrame("stokesI"), zerow=False)
    print(block_vis)
    vis = convert_blockvisibility_to_visibility(block_vis)
    
    vis.data['uvw'][..., 2] = 0.0
    
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
    
    advice = advise_wide_field(vis, guard_band_image=1.0, delA=0.02)
    
    pb_npixel = 2048
    d2r = numpy.pi / 180.0
    pb_cellsize = d2r * FOV_deg / pb_npixel
    
    cellsize = advice['cellsize']
    npixel = 512
    
    if show:
        plt.clf()
        plt.plot(-vis.u, -vis.v, '.', color='b', markersize=0.2)
        plt.plot(vis.u, vis.v, '.', color='b', markersize=0.2)
        plt.xlabel('U (wavelengths)')
        plt.ylabel('V (wavelengths)')
        plt.title('UV coverage')
        plt.savefig('uvcoverage.png')
        plt.show(block=False)
        
        ha = numpy.pi * block_vis.time / 43200.0
        dec = phasecentre.dec.rad
        latitude = block_vis.configuration.location.lat.rad
        az, el = hadec_to_azel(ha, dec, latitude)
        plt.clf()
        plt.plot(ha, az, '.', label='Azimuth (rad)')
        plt.plot(ha, el, '.', label='Elevation (rad)')
        plt.xlabel('HA (rad)')
        plt.legend()
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
        from rascil.processing_components.simulation.testing_support import create_test_skycomponents_from_s3
        
        original_components = create_test_skycomponents_from_s3(flux_limit=flux_limit,
                                                                phasecentre=phasecentre,
                                                                polarisation_frame=PolarisationFrame("stokesI"),
                                                                frequency=numpy.array(frequency),
                                                                radius=pb_cellsize * pb_npixel / 2.0)
        # Primary beam points to the phasecentre
        offset_direction = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')
    
    # Uniform weighting
    psf = create_image_from_visibility(vis, npixel=npixel, frequency=frequency,
                                       nchan=nfreqwin, cellsize=cellsize, phasecentre=phasecentre,
                                       polarisation_frame=PolarisationFrame("stokesI"))

    if use_natural:
        print("Using natural weighting")
    else:
        print("Using uniform weighting")
        vis = weight_list_serial_workflow([vis], [psf])[0]
        block_vis = convert_visibility_to_blockvisibility(vis)
    
    print("Inverting to get on-source PSF")
    psf_list = invert_list_rsexecute_workflow([vis], [psf], '2d', dopsf=True)
    psf, sumwt = rsexecute.compute(psf_list, sync=True)[0]
    export_image_to_fits(psf, 'PSF_rascil.fits')
    if show:
        show_image(psf, cm='gray_r', title='PSF', vmin=-0.01, vmax=0.1)
        plt.savefig('PSF_rascil.png')
        plt.show(block=False)
    
    model = create_image_from_visibility(vis, npixel=npixel, frequency=frequency,
                                         nchan=nfreqwin, cellsize=cellsize, phasecentre=offset_direction,
                                         polarisation_frame=PolarisationFrame("stokesI"))
    
    # ### Calculate the voltage pattern without pointing errors
    vp = create_image_from_visibility(block_vis, npixel=pb_npixel, frequency=frequency,
                                      nchan=nfreqwin, cellsize=pb_cellsize, phasecentre=phasecentre,
                                      override_cellsize=False)
    
    if show:
        pb = create_pb(vp, pbtype, pointingcentre=phasecentre, use_local=not use_radec)
        print("Primary beam:", pb)
        show_image(pb, title='%s: primary beam' % context)
        plt.savefig('PB_rascil.png')
        export_image_to_fits(pb, 'PB_rascil.fits')
        plt.show(block=False)
    
    print("Constructing voltage pattern")
    vp = create_vp(vp, pbtype, pointingcentre=phasecentre, use_local=not use_radec)
    print("Voltage pattern:", vp)
    pt = create_pointingtable_from_blockvisibility(block_vis)
    
    no_error_pt = simulate_pointingtable(pt, 0.0, 0.0)

    export_pointingtable_to_hdf5(no_error_pt, 'pointingsim_%s_noerror_pointingtable.hdf5' % context)
    no_error_gt = simulate_gaintable_from_pointingtable(block_vis, original_components, no_error_pt, vp,
                                                        use_radec=use_radec)
    # Each component in original components becomes a separate skymodel
    no_error_sm = [SkyModel(components=[original_components[i]], gaintable=no_error_gt[i])
                   for i, _ in enumerate(original_components)]
    
    # We do this in chunks of eight to avoid creating all visibilities at once
    no_error_blockvis = copy_visibility(block_vis, zero=True)
    
    print("Predicting error-free visibilities in chunks of %d skymodels" % ngroup)
    future_vis = rsexecute.scatter(no_error_blockvis)
    chunks = [no_error_sm[i:i + ngroup] for i in range(0, len(no_error_sm), ngroup)]
    for chunk in chunks:
        temp_vis = predict_skymodel_list_compsonly_rsexecute_workflow(future_vis, chunk, context='2d', docal=True)
        work_vis = rsexecute.compute(temp_vis, sync=True)
        for w in work_vis:
            no_error_blockvis.data['vis'] += w.data['vis']
        assert numpy.max(numpy.abs(no_error_blockvis.data['vis'])) > 0.0
    
    no_error_vis = convert_blockvisibility_to_visibility(no_error_blockvis)
    print("Inverting to get dirty image")
    dirty_list = invert_list_rsexecute_workflow([no_error_vis], [model], '2d')
    dirty, sumwt = rsexecute.compute(dirty_list, sync=True)[0]
    print(qa_image(dirty))
    export_image_to_fits(dirty, 'dirty_rascil.fits')
    if show:
        show_image(dirty, cm='gray_r', title='Dirty image')  # , vmin=-0.01, vmax=0.1)
        plt.savefig('dirty_rascil.png')
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
    for pe in pes:
    
        numpy.random.seed(seed)

        result = dict()
        result['context'] = context
        result['nb_name'] = sys.argv[0]
        result['plotfile'] = plotfile
        result['hostname'] = socket.gethostname()
        result['epoch'] = epoch
        basename = os.path.basename(os.getcwd())
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

        if time_series is '':
            error_pt = simulate_pointingtable(pt, pointing_error=pointing_error * a2r,
                                              static_pointing_error=static_pointing_error * a2r,
                                              global_pointing_error=global_pointing_error * a2r)
        else:
            error_pt = simulate_pointingtable_from_timeseries(pt, type=time_series)

        export_pointingtable_to_hdf5(error_pt,
                                     'pointingsim_%s_error_%.0farcsec_pointingtable.hdf5' % (context, pe))
        
        print("Scaling pointing errors in %.3f, %.3f" % (scale[0], scale[1]))
        error_pt.pointing[..., 0] *= scale[0]
        error_pt.pointing[..., 1] *= scale[1]
        
        error_gt = simulate_gaintable_from_pointingtable(block_vis, original_components, error_pt, vp,
                                                         use_radec=use_radec)
        
        error_sm = [SkyModel(components=[original_components[i]], gaintable=error_gt[i])
                    for i, _ in enumerate(original_components)]
        
        error_blockvis = copy_visibility(block_vis, zero=True)
        
        print("Predicting corrupted visibilities in chunks of %d skymodels" % ngroup)
        future_vis = rsexecute.scatter(error_blockvis)
        chunks = [error_sm[i:i + ngroup] for i in range(0, len(error_sm), ngroup)]
        for chunk in chunks:
            temp_vis = predict_skymodel_list_compsonly_rsexecute_workflow(future_vis, chunk, context='2d', docal=True)
            work_vis = rsexecute.compute(temp_vis, sync=True)
            for w in work_vis:
                error_blockvis.data['vis'] += w.data['vis']
            assert numpy.max(numpy.abs(error_blockvis.data['vis'])) > 0.0
        error_blockvis.data['vis'] -= no_error_blockvis.data['vis']
        
        error_blockvis = addnoise_visibility(error_blockvis, tsys)
        
        print("Inverting to get on-source dirty image")
        error_vis = convert_blockvisibility_to_visibility(error_blockvis)
        dirty_list = invert_list_rsexecute_workflow([error_vis], [model], '2d')
        dirty, sumwt = rsexecute.compute(dirty_list, sync=True)[0]
        export_image_to_fits(dirty, 'PE_%.1f_arcsec_rascil.fits' % pe)
        if show:
            show_image(dirty, cm='gray_r', title='Pointing error %.1f arcsec rascil' % pe)
            plt.savefig('PE_%.1f_arcsec_rascil.png' % pe)
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
            outlier_dirty = invert_list_rsexecute_workflow([error_vis], [outlier_model], '2d')
            outlier_dirty, outlier_sumwt = rsexecute.compute(outlier_dirty, sync=True)[0]
            
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
    for ifield, field in enumerate(['onsource_maxabs', 'onsource_rms', 'onsource_medianabs']):
        plt.loglog(pes, [result[field] for result in results], '-', label=field, color=colors[ifield])
    
    if outlier:
        for ifield, field in enumerate(['outlier_maxabs', 'outlier_rms', 'outlier_medianabs']):
            plt.loglog(pes, [result[field] for result in results], '--', label=field, color=colors[ifield])
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
