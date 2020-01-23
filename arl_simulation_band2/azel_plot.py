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
import os
import sys
import time

from rascil.data_models.parameters import rascil_path

results_dir = rascil_path('test_results')

import numpy

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from rascil.data_models.polarisation import PolarisationFrame

from rascil.processing_components.visibility.base import create_blockvisibility
from rascil.processing_components.simulation.configurations import create_configuration_from_MIDfile
from processing_library.util.coordinate_support import hadec_to_azel

import logging

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

import pprint

pp = pprint.PrettyPrinter()

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
                        help='s3sky or singlesource or null')
    
    parser.add_argument('--rmax', type=float, default=1e5,
                        help='Maximum distance of station from centre (m)')
    
    parser.add_argument('--global_pe', type=float, nargs=2, default=[0.0, 0.0], help='Global pointing error')
    parser.add_argument('--static_pe', type=float, nargs=2, default=[0.0, 0.0], help='Multipliers for static errors')
    parser.add_argument('--dynamic_pe', type=float, default=1.0, help='Multiplier for dynamic errors')
    parser.add_argument('--nnodes', type=int, default=1, help='Number of nodes')
    parser.add_argument('--nthreads', type=int, default=1, help='Number of threads')
    parser.add_argument('--memory', type=int, default=8, help='Memory per worker (GB)')
    parser.add_argument('--nworkers', type=int, default=8, help='Number of workers')
    parser.add_argument('--flux_limit', type=float, default=1.0, help='Flux limit (Jy)')
    parser.add_argument('--show', type=str, default='False', help='Show images?')
    parser.add_argument('--export_images', type=str, default='False', help='Export images in fits format?')
    parser.add_argument('--ngroup_visibility', type=int, default=8, help='Process in visibility groups this large')
    parser.add_argument('--ngroup_components', type=int, default=8, help='Process in component groups this large')
    parser.add_argument('--npixel', type=int, default=512, help='Number of pixels in image')
    parser.add_argument('--seed', type=int, default=18051955, help='Random number seed')
    parser.add_argument('--snapshot', type=str, default='False', help='Do snapshot only?')
    parser.add_argument('--opposite', type=str, default='False',
                        help='Move source to opposite side of pointing centre')
    parser.add_argument('--offset_dir', type=float, nargs=2, default=[0.0, 0.0], help='Multipliers for null offset')
    parser.add_argument('--pbradius', type=float, default=2.0, help='Radius of sources to include (in HWHM)')
    parser.add_argument('--pbtype', type=str, default='MID', help='Primary beam model: MID or MID_GAUSS')
    parser.add_argument('--use_agg', type=str, default="True", help='Use Agg matplotlib backend?')
    parser.add_argument('--declination', type=float, default=-45.0, help='Declination (degrees)')
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
    
    declination = args.declination
    use_radec = args.use_radec == "True"
    use_natural = args.use_natural == "True"
    export_images = args.export_images == "True"
    time_series = args.time_series
    pointing_file = args.pointing_file
    scale = numpy.array(args.scale)
    tsys = args.tsys
    integration_time = args.integration_time
    time_range = args.time_range
    time_chunk = args.time_chunk
    snapshot = args.snapshot == 'True'
    opposite = args.opposite == 'True'
    offset_dir = args.offset_dir
    pbtype = args.pbtype
    pbradius = args.pbradius
    reference_pointing = args.reference_pointing == "True"
    pointing_directory = args.pointing_directory
    rmax = args.rmax
    flux_limit = args.flux_limit
    npixel = args.npixel
    shared_directory = args.shared_directory
    
    seed = args.seed
    print("Random number seed is", seed)
    show = args.show == 'True'
    context = args.context
    nworkers = args.nworkers
    nnodes = args.nnodes
    threads_per_worker = args.nthreads
    memory = args.memory
    ngroup_visibility = args.ngroup_visibility
    ngroup_components = args.ngroup_components
    
    basename = os.path.basename(os.getcwd())
    
    time_started = time.time()
    
    # Set up details of simulated observation
    nfreqwin = 1
    diameter = 15.0
    frequency = [1.36e9]
    channel_bandwidth = [1e7]
    
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=declination * u.deg, frame='icrs', equinox='J2000')
    location = EarthLocation(lon="21.443803", lat="-30.712925", height=0.0)
    
    # Do each 30 minutes in parallel
    start_times = numpy.arange(time_range[0] * 3600, time_range[1] * 3600, time_chunk)
    end_times = start_times + time_chunk
    print("Start times for chunks:")
    pp.pprint(start_times)
    
    
    def valid_elevation(time, location, phasecentre):
        ha = numpy.pi * time / 43200.0
        dec = phasecentre.dec.rad
        az, el = hadec_to_azel(ha, dec, location.lat.rad)
        return el > 15.0 * numpy.pi / 180.0
    
    
    number_valid_times = 0
    valid_start_times = []
    for t in start_times:
        if valid_elevation(t, location, phasecentre) or valid_elevation(t + time_chunk, location, phasecentre):
            valid_start_times.append(t)
            number_valid_times += 1
    
    assert number_valid_times > 0, "No data above elevation limit"
    
    print("Start times for chunks above elevation limit:")
    start_times = valid_start_times
    pp.pprint(start_times)
    
    times = [numpy.arange(start_times[itime], end_times[itime], integration_time) for itime in
             range(len(start_times))]
    print("Observation times:")
    pp.pprint(times)
    
    s2r = numpy.pi / (12.0 * 3600)
    rtimes = s2r * numpy.array(times)
    ntimes = len(rtimes.flat)
    nchunks = len(start_times)
    
    assert ntimes > 0, "No data above elevation limit"
    
    print('%d integrations of duration %.1f s processed in %d chunks' % (ntimes, integration_time, nchunks))
    
    mid = create_configuration_from_MIDfile('%s/ska1mid_local.cfg' % shared_directory, rmax=rmax,
                                            location=location)
    
    bvis_list = [create_blockvisibility(mid, rtimes[itime], frequency=frequency,
                                        channel_bandwidth=channel_bandwidth, weight=1.0,
                                        phasecentre=phasecentre,
                                        polarisation_frame=PolarisationFrame("stokesI"),
                                        zerow=True)
                 for itime in range(nchunks)]
    
    plt.clf()
    r2d = 180.0 / numpy.pi
    for ivis in range(nchunks):
        bvis = bvis_list[ivis]
        ha = numpy.pi * bvis.time / 43200.0
        dec = phasecentre.dec.rad
        latitude = bvis.configuration.location.lat.rad
        az, el = hadec_to_azel(ha, dec, latitude)
        if ivis == 0:
            plt.plot(bvis.time, r2d * az, '.', color='r', label='Azimuth (deg)')
            plt.plot(bvis.time, r2d * el, '.', color='b', label='Elevation (deg)')
        else:
            plt.plot(bvis.time, r2d * az, '.', color='r')
            plt.plot(bvis.time, r2d * el, '.', color='b')
    plt.xlabel('HA (s)')
    plt.ylabel('Angle')
    plt.legend()
    plt.title('Azimuth and elevation vs hour angle')
    plt.savefig('azel.png')
    plt.show(block=False)
