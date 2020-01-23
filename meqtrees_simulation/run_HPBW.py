#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This shows how to run a TDL script in a pipeline (aka batch-mode, aka headless mode)
#
import sys, numpy as np
from Timba.Apps import meqserver
from Timba.TDL import Compile
from Timba.TDL import TDLOptions
import os
import numpy as np
import h5py

import h5py

def runit(pe, dirname):
    # This starts a meqserver. Note how we pass the "-mt 2" option to run two threads.
    # A proper pipeline script may want to get the value of "-mt" from its own arguments (sys.argv).
    print "Starting meqserver";
    mqs = meqserver.default_mqs(wait_init=10,extra=["-mt","16"]);

    # Once we're connected to a server, some cleanup is required before we can exit the script.
    # Since we want to perform this cleanup regardless of whether the script ran to completion
    # or was interrupted by an exception midway through, we use a try...finally block.
    try:

        if pe==0.0:
            print 'im 0.0'
            lerrs = np.array([0.0])
            merrs = np.array([0.0])
            ms_name = dirname+'PE_0.0_arcsec.ms'
            os.system('rm -rf '+dirname+ms_name)
            os.system('cp -r '+dirname+'setupvis.ms '+ms_name)

        else: 
            ms_name = dirname+'PE_%s_arcsec.ms'%str(pe)
            os.system('rm -rf '+ms_name)
            os.system('cp -r '+dirname+'setupvis.ms '+ms_name)
            lerrs, merrs = hdf2npy(pe)
          #  lerrs = np.load(dirname+'RA_offsets%s.npy'%pe)
          #  merrs = np.load(dirname+'Dec_offsets%s.npy'%pe)
          #  lerrs*=206265
          #  merrs*=206265
           # print lerrs, merrs
        #np.save('/usr/local/lib/python2.7/dist-packages/Cattery/save_sigma', float(pe))
        source_list_name = dirname+'source_list.lsm.html'
        print 'making ', ms_name
        print '========== Making config file';
        make_config_file(lerrs,merrs,ms_name, source_list_name);
        TDLOptions.config.read("PEs.tdl.conf");
        script = "turbo-sim.py";
        print "========== Compiling batch job 1";
        mod,ns,msg = Compile.compile_file(mqs,script, config = "turbo-sim");
        print "========== Running batch job 1";
        mod._simulate_MS(mqs,None,wait=True);
    except: 
        print 'Corrupted MS not made properly'
        
       
  ### Cleanup time
    finally:
        print "Stopping meqserver";
     # this halts the meqserver
        meqserver.stop_default_mqs();
    # now we can exit
        print "Bye!";


def hdf2npy(pe):
    if 1:
        file_name = 'pointingsim_error_arcsec%s_pointingtable.hdf5'%str(pe)
        print 'loading PEs from ', file_name
        f = h5py.File(file_name)
        group = f['PointingTable0']
        d = group['data'].value
        d2 = np.array(d)
        for j in range(len(d)):
            try:
                dnpy = np.dstack((dnpy,d[j][0][:,0,0,:]))
            except:
                dnpy = d[j][0][:,0,0,:]

       

        for i in range(65):
            try:
                allerrs = np.vstack((allerrs,dnpy[:,:,i]))
            except:
                allerrs = dnpy[:,:,i]
        print allerrs.shape
        lerrs = allerrs[:,0].astype(np.float)                    ##### rascil PEs, in radians
        merrs= allerrs[:,1].astype(np.float)
        lerrs*=206265
        merrs*=206265
        print lerrs.shape, merrs.shape
        f.close()
        return lerrs, merrs




def make_config_file(lerrs, merrs,ms_name, source_list_name):
    if len(lerrs) ==1:
        print 'zero ', lerrs
        lerrs = 0
        merrs = 0
    else:
        lerrs =     ' '.join(map(str,lerrs))
        merrs =     ' '.join(map(str,merrs))
    print 'l and m ', lerrs, merrs
    os.system('rm -rf PEs.tdl.conf')
    f = open("PEs.tdl.conf","w")
    f.write('[turbo-sim]\n')
    f.write('# compile-time options follow\n')
    f.write('ms_sel.msname = %s\n'%ms_name)
    f.write('ms_sel.ms_ifr_subset_str = all\n')
    f.write('ms_sel.ms_corr_sel = 1\n')
    f.write('run_purr = 0\n')
    f.write('sim_mode = sim only\n')
    f.write('read_ms_model = 0\n')
    f.write('tensormeqmaker.psv_class = PSVTensor\n')
    f.write('me.use_jones_inspectors = 1\n')
    f.write('me.use_skyjones_visualizers = 0\n')
    f.write('me.use_smearing = 0\n')
    f.write('uvw_source = from MS\n')
    f.write('uvw_refant = default\n')
    f.write('me.sky.tiggerskymodel = 1\n')
    f.write('tiggerlsm.filename = %s\n'%source_list_name)
    f.write('tiggerlsm.lsm_subset = all\n')
    f.write('tiggerlsm.null_subset = None\n')
    f.write('tiggerlsm.solvable_sources = 0\n')
    f.write('me.sky.siamese_oms_gridded_sky = 0\n')
    f.write('me.sky.siamese_agw_azel_sky = 0\n')
    f.write('me.sky.siamese_oms_transient_sky = 0\n')
    f.write('me.sky.siamese_oms_fitsimage_sky = 0\n')
    f.write('me.ncorr_enable = 0\n')
    f.write('me.z_enable = 0\n')
    f.write('me.l_enable = 0\n')
    f.write('me.e_enable = 1\n')
    f.write('me.e_module = Siamese_OMS_analytic_beams\n')
    f.write('model_type = use_circular_aperture\n')
    f.write('analytic_beams.circular_aperture_beam.bf = 1.0\n')
    f.write('analytic_beams.circular_aperture_beam.dish_sizes = 15\n')

    f.write('me.epe_enable = 1\n')
    f.write('oms_pointing_errors.station_subset = all\n')
    f.write('oms_pointing_errors.pe_l.error_model = ListOfValues\n')
    f.write('oms_pointing_errors.pe_l.values_str = %s\n'%lerrs)
    f.write('oms_pointing_errors.pe_m.error_model = ListOfValues\n')
    f.write('oms_pointing_errors.pe_m.values_str = %s\n'%merrs)
    f.write('me.e_advanced = 0\n')
    f.write('me.p_enable = 0\n')
    f.write('me.d_enable = 0\n')
    f.write('me.g_enable = 0\n')
    f.write('me.ip_enable = 0\n')
    f.write('noise_stddev = None\n')
    f.write('noise_from_sefd = 0\n')
    f.write('random_seed = time\n')
    f.write('# runtime options follow\n')
    f.write('img_sel.image_clean_gain = 0.1\n')
    f.write('img_sel.image_clean_method = clark\n')
    f.write('img_sel.image_clean_niter = 1000\n')
    f.write('img_sel.image_clean_resetmodel = 1\n')
    f.write('img_sel.image_clean_threshold = 0Jy\n')
    f.write('img_sel.image_viewer = tigger\n')
    f.write('img_sel.imaging_arcmin = 30.0\n')
    f.write('img_sel.imaging_chanmode = 1 (average all)\n')
    f.write('img_sel.imaging_column = MODEL_DATA\n')
    f.write('img_sel.imaging_custom_ms_select = 0\n')
    f.write('img_sel.imaging_enable_wproj = 0\n')
    f.write('img_sel.imaging_freqmode = frequency\n')
    f.write('img_sel.imaging_ifrs = all\n')
    f.write('img_sel.imaging_npix = 18703\n')
    f.write('img_sel.imaging_padding = 10.0\n')
    f.write('img_sel.imaging_phasecenter = default\n')
    f.write('img_sel.imaging_stokes = I\n')
    f.write('img_sel.imaging_taper_gauss = 0\n')
    f.write('img_sel.imaging_weight = uniform\n')
    f.write('img_sel.output_fitsname = default\n')
    f.write('ms_sel.ms_taql_str = None\n')
    f.write('ms_sel.output_column = MODEL_DATA\n')
    f.write('ms_sel.select_channels = 0\n')
    f.write('ms_sel.tile_size = 65\n')
    f.close()




  




if __name__ == '__main__':
    import sys, numpy as np
    from Timba.Apps import meqserver
    from Timba.TDL import Compile
    from Timba.TDL import TDLOptions
     
    #errs = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0,64.0, 128.0, 256.0]
    #for i in errs:
    i = np.float(sys.argv[1])
    print i
    dirname= sys.argv[2]
    runit(i, dirname)









