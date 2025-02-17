# Author P Hartley March 2019
# 
# notes: had to change imclean in simutil.py so that minbp = 0. and pbcor = pbcor are passed to clean
# rms 30 arcsec for pointing offsets
# config files: alma_cycle1_1.cfg
#               ska1mid_modified_no_meerkat.cfg - homegeneous, and keeps ALMA PB for now
# SKA-mid resolution: 

from matplotlib import pyplot as plt
import sys
sys.path.append("/home/p.hartley/casa-release-5.4.1-32.el7/analysis_scripts")# move this
import analysisUtils as au
from astropy.io import fits
import numpy as np


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



config_file = '/home/p.hartley/Documents/CASA_sims/casa_scripts/ska1mid.cfg' #treat as homogeneous
print  '.'+config_file.split('/')[-1].strip('cfg')
config_file_prefix = '.'+config_file.split('/')[-1].strip('cfg')
print config_file
myVPim = '/home/p.hartley/Documents/CASA_sims/beam_conversion/convert-to-casa-images/results/main_beam_x-directivity.image'
mycomplexvp = '/home/p.hartley/Documents/CASA_sims/beam_conversion/convert-to-casa-images/results_vp/main_beam_FrEnd_uv-voltage_pattern.image'
#mycomplexvp = '/home/p.hartley/Documents/CASA_sims/beam_conversion/convert-to-casa-images/results_bak/main_beam_FrEnd-el_field_pattern.image'




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def source_list(nsources,gap, centre, s, dirname, layout = 'grid'):
    os.system('rm -rf points_grid.cl')
    cl.done()
    f = open(dirname+"source_list.txt","w")
    f.write('#format: name ra_h ra_m ra_s dec_d dec_m dec_s i\n') # write coords to ascii for meqtrees use
    count = 0
    lamb = 3e8/(1.4e9)
    HPBW = (70*lamb)/15.0 #in degrees
      #  HPBW = HPBW/360.*2*np.pi
    HPBW/=2.  
    print 'HPBW halved: ',HPBW
    HPBW = (HPBW/360.)*2*np.pi
    if layout =='grid':
        gridpoints = np.arange(nsources)-(nsources/2)
        print gridpoints
        for i in gridpoints:
            for j in gridpoints:
                ra_shift =  i*gap/rad2arcsec
                dec_shift = j*gap/rad2arcsec
                ra = centre['m0']['value']+ra_shift
                dec = centre['m1']['value']+dec_shift # need to change to absolute distance
                radec = au.rad2radec(ra, dec, hmsdms = True)
                f.write('{} {} {} {} {} {} {} {}\n'.format(count,radec.split(', ')[0].split('h')[0],\
                           radec.split(', ')[0].split('h')[1].split('m')[0], \
                           radec.split(', ')[0].split('m')[1].split('s')[0],\
                           radec.split(', ')[1].split('d')[0],\
                           radec.split(', ')[1].split('d')[1].split('m')[0], \
                           radec.split(', ')[1].split('m')[1].split('s')[0], s ))
                cl.addcomponent(dir="J2000 {} {}".format(radec.split(', ')[0],radec.split(', ')[1]), flux=s, fluxunit='Jy', freq='1.4GHz', shape="point") # make CASA comp list
                count+=1
    if layout == 'single': # ! bug in simobserve that won't make image with only 1 comp - have fixed but need to thread a switch through
        ra = centre['m0']['value']
        dec = centre['m1']['value']+HPBW
        radec = au.rad2radec(ra, dec, hmsdms = True)    
        f.write('{} {} {} {} {} {} {} {}\n'.format(count,radec.split(', ')[0].split('h')[0],\
                           radec.split(', ')[0].split('h')[1].split('m')[0], \
                           radec.split(', ')[0].split('m')[1].split('s')[0],\
                           radec.split(', ')[1].split('d')[0],\
                           radec.split(', ')[1].split('d')[1].split('m')[0], \
                           radec.split(', ')[1].split('m')[1].split('s')[0], s ))
        cl.addcomponent(dir="J2000 {} {}".format(radec.split(', ')[0],radec.split(', ')[1]), flux=s, fluxunit='Jy', freq='1.4GHz', shape="point") # make CASA comp list
    cl.rename(dirname+'points_grid.cl')
    cl.close()
    cl.done()
            

def do_plot_azel(azs, els):
        plt.clf()
        plt.scatter(azs,els)
        plt.savefig('az_el_plot.png')
     

    

def observe_simobserve(projectname, config_file, layout):
    ms_direction = 'J2000 23h00m00.031s -40d59m59.6s'
    d =  me.direction(ms_direction.split()[0], ms_direction.split()[1], ms_direction.split()[2])# prepare to convert to degrees f
    if make_source_list:
        layout = 'single'
        if layout =='single':
            make_modsize = True
        centre = d
        source_list(3, 5*60, centre, 1.0, dirname, layout = 'grid')
    # generate some idealised visibilities using the simobserve task
    simobserve(project=projectname, 
                        complist = 'points_grid.cl',
                        compwidth='128MHz' , 
                        setpointings=True, 
                        integration='600s',  
                        direction=ms_direction,  
                        mapsize= '3arcsec', # mapsize < PB forces single pointing
                        obsmode='int', 
                        antennalist='../ska1mid_modified.cfg', #config_file, 
                        hourangle='transit', 
                        totaltime='600s',  
                        thermalnoise='', 
                        graphics='both', 
                        overwrite=True,
                        )
                        
def observe_simulator(projectname, config_file, NUM_SPW, fieldDir, dirname, layout) :
    make_source_list = 1
    if make_source_list:
        if layout =='single':
            make_modsize = True
        centre = fieldDir
        source_list(5, 15*60, centre, 1.0, dirname , layout = layout)
    DUMP_TIME = 1.                 # integration time in seconds
    SCAN_TIME = (DUMP_TIME*1.)/3600.    # make each scan last a few integrations, units of hours
    NUM_SCANS = 65             # accumulate this many scans (use odd number)
    HA_RANGE = 8.0                      # total range of HA coverage in hours (use > 0.0)
    FREQ = 1.4                          # centre frequency in GHz
    FREQ_FBW = 0.01                     # frequency fractional bandwidth per channel
    NUM_CHAN = 1                      # number of spectral channels
                       # number of spectral windows (use odd number)
    FREQ_SPW = 0.014 
    REF_TIME = me.epoch('UTC','50000.0d')
    xfreq = str(FREQ) + 'GHz'
    FREQ_RES = str(FREQ*FREQ_FBW*1000.) +'MHz' 
          
    tel_name = 'SKA1-MID'
    output_ms = projectname+'.'+config_file.split('/')[-1].strip('cfg')+'ms'
    print 'output_ms: ', output_ms
    os.system( 'rm -rf ' + output_ms )
    output_table = projectname+'.'+config_file.split('/')[-1].strip('cfg')+'-vp.tab'
    os.system( 'rm -rf ' + output_table)
    sm.open( output_ms )
  # read configuration file
    print 'open'
    f = open( config_file, 'r' )
    xx = []; yy = []; zz = []; dd = []; ss = []
    while True:
      line = f.readline()
      if not line: break
      items = line.split()
      if ( items[0] != '#' ):
        xx.append( float( items[0] ) )
        yy.append( float( items[1] ) )
        zz.append( float( items[2] ) )
        dd.append( float( 15.0 ) )
        ss.append( str( items[4] ) )
    f.close()
    vp.reset()
    vp.setpbairy(telescope=tel_name,dishdiam='15m',maxrad='5deg')
   # print 'setting vp from ', mycomplexvp
   # vp.setpbimage(telescope='SKA1-MID', compleximage=mycomplexvp,antnames='*')
    print 'saving vp table in: ', output_table
    vp.saveastable(tablename = output_table)
    #posskamid = me.observatory( 'SKA1-MID')
    posskamid = me.position('wgs84', '21.443803deg', '-30.712925deg', '1053.000000m')
    sm.setconfig(telescopename = tel_name, x = xx, y = yy, z = zz, dishdiameter = dd, antname = ss, mount = 'ALT-AZ',  coordsystem = 'global',referencelocation = posskamid  )
    jvar = NUM_SPW +1
    for j in range(1, jvar):
    # define the spectral windows
        ifnum = (j - 1)
        xspwname = 'IF' + str(ifnum)
        spwfreq = str(FREQ + float( j - ((NUM_SPW + 1) / 2) )*FREQ_SPW) + 'GHz'
        sm.setspwindow( spwname = xspwname, freq = spwfreq, deltafreq = FREQ_RES, freqresolution = FREQ_RES, nchannels = NUM_CHAN, stokes = 'XX' ) 
        print (xspwname, spwfreq)
    print 'config set'
    sm.setfeed( mode = 'perfect X Y', pol = [''] )
    sm.setfield( sourcename = 'ptcentre', sourcedirection = fieldDir )
    sm.setlimits( shadowlimit = 0.001, elevationlimit = '15.0deg' )
    sm.setauto( autocorrwt = 0.0 )
    sm.settimes( integrationtime = str( DUMP_TIME ) + 's', usehourangle = True, referencetime = REF_TIME )       
    ivar = NUM_SCANS +1
    for i in xrange(1, ivar):
    # set the start and end times.
        if (NUM_SCANS > 1): 
            begin = float( i - ((NUM_SCANS + 1) / 2.) )/float(NUM_SCANS-1)*HA_RANGE  -  SCAN_TIME / 2.
        else:
            begin =  - SCAN_TIME / 2.    
        end = begin +  SCAN_TIME 
        # convert to string.
        begin = str( begin ) + 'h'
        end = str( end ) + 'h'
        for j in xrange(1, jvar):
            ifnum = (j - 1)
            xspwname = 'IF' + str(ifnum)
            sm.observe( sourcename= 'ptcentre', spwname= xspwname, starttime = begin, stoptime = end )
  # get model image from fits file 
    sm.setvp( dovp = True, usedefaultvp = False, vptable = output_table, dosquint = False)
  # generate the predicted visibilities.
    sm.predict( complist = dirname+'points_grid.cl')
    print 'predicted'
    sm.corrupt()
    print 'corrupted'

    sm.close()
   
             
             
                       
                      
                        
def corrupt(projectname, config_file_prefix,offset, dirname,PEs = False,  plot_azel=False):   
    if PEs:
        offset_name = np.copy(offset)
        print 'in PEs'
        #### sm.setpointingerror() 
        # no documentation for how to contruct a PE E Jones table and apply it:
        # turns out that the setpointingerror function is not written 
        # instead using the pointings from the POINTINGS table and modifying them manually
    
        
        # overwrite the idealised visibilities with corrupted data
        # pointing and antenna corruptions are applied before the visibilities are obtained
        # other corruptions are applied to the visibilities themselves                     
        tb.open(projectname+config_file_prefix+'ms/ANTENNA')
        pos=tb.getcol('POSITION')
        ant=tb.getcol('NAME')
        pad=tb.getcol('STATION')
        tb.close()
        
        # open pointings file to modify    
        tb.open(projectname+config_file_prefix+'ms/POINTING', nomodify = False)
        data = tb.getcol('DIRECTION')
        antenna = tb.getcol('ANTENNA_ID')
        time = tb.getcol('TIME')
        dataRA = data[0,0,:]
        dataDec = data[1,0,:]
        azimuths  = np.array([])
        elevations = np.array([])

        # apply a random static offset to each antenna in polar coords and decompose to az el
        # example error: rms 30 arcsec = 0.5arcmin not 0.5deg!
   #     ant_offsets = np.zeros(len(ant))+0.5/60.# offsets in degrees here for plotting
        ant_offsets =np.zeros(len(ant))+(0.5/60.)#np.random.normal(0, 0./60., len(ant)) # offsets in degrees here for plotting
        print 'len dataRA = ', len(dataRA)
        time_offsets = np.zeros(len(dataRA))#np.random.normal(0,60/3600., len(dataRA))
        np.random.seed(0)
        radec_offsets = np.random.normal(0,1, len(dataRA))
        print 'offset rands: ', radec_offsets
        radec_offsets*= (offset/3600.)
        
        #ant_thetas = ((2*np.pi)/len(ant))*(np.arange(len(ant))+1)##np.random.rand(len(ant))*2*np.pi#2*np.pi)/(np.arange(len(ant))+1)#  
        ant_thetas = np.zeros(len(ant))#)))np.random.rand(len(ant))*2*np.pi
        time_thetas = np.zeros(len(dataRA))#)p.random.rand(len(dataRA))*2*np.pi
 #np    .zeros(len(ant))#
        np.random.seed(2)
        radec_thetas = np.random.rand(len(dataRA))
        print 'theta rands: ', radec_thetas
        radec_thetas*= 2*np.pi
        offset/=3600.
        np.random.seed(5)
        ra_offsets = ((np.random.normal(0,1,len(dataRA)))*offset)
        np.random.seed(6)
        dec_offsets = (np.random.normal(0,1,len(dataRA)))*offset
        ra_offsets = (ra_offsets/360.)*2*pi
        dec_offsets = (dec_offsets/360.)*2*pi
        print dec_offsets.shape
        if 0:
            ant_az_offsets = np.cos(ant_thetas)*ant_offsets 
            ant_el_offsets = np.sin(ant_thetas)*ant_offsets
            time_az_offsets = np.cos(time_thetas)*time_offsets 
            time_el_offsets = np.sin(time_thetas)*time_offsets
            ra_offsets = (np.cos(radec_thetas))*radec_offsets
            dec_offsets = (np.sin(radec_thetas))*radec_offsets
        print ra_offsets
        print dec_offsets
        print 'ra dec ave', np.mean(ra_offsets), np.mean(dec_offsets)   
        print 'ra dec median', np.median(ra_offsets), np.median(dec_offsets)    
        pe32 = np.load('arcsec%s.npy'%offset_name)
        print pe32.shape
        for i in range(65):
            try:
                pearr = np.vstack((pearr,pe32[:,:,i]))
            except:
                pearr = pe32[:,:,i]
        ra_offsets = pearr[:,0].astype(np.float)                    ##### rascil PEs
        dec_offsets= pearr[:,1].astype(np.float)

     #   plt.clf()
     #   print ant_az_offsets.shape
     #   print ant_el_offsets.shape
     #   plt.scatter(ant_az_offsets, ant_el_offsets, color = '#ac0c4b', linewidth = 0.0)
     #   plt.axhline(0, color = 'black')
     #   plt.axvline(0, color = 'black')
     #   plt.scatter(np.mean(ant_az_offsets),np.mean(ant_el_offsets) )   
     #   plt.ylabel('elevation /degrees')
     #   plt.xlabel('azimuth /degrees')
     #   plt.savefig(dirname+'pointing_offsets.png') 

        
        az_tot = 0
        el_tot = 0
        RA_tot = 0
        Dec_tot = 0
        RA_offsets = np.array([])
        Dec_offsets = np.array([])
        # use measures to convert ra dec to az el
        print 'number of visibilities: ', len(dataRA)
        for scan in xrange(len(dataRA)):# note len(dataRA) not equal to n visibilities if dump<scan time
            a = me.epoch('UTC',str(time[scan])+'s')      # a time  
            me.measure(a, 'tai') # convert to IAT    
            me.doframe(a) # set time in frame   
            antpos=me.position('ITRF',qa.quantity(pos[0,antenna[scan]],'m'),qa.quantity(pos[1,antenna[scan]],'m')\
                ,qa.quantity(pos[2,antenna[scan]],'m'))#check ITRF
            me.doframe(antpos)
            b =  me.direction('j2000', str(dataRA[scan])+'rad', str(dataDec[scan])+'rad' )
            az_el =  me.measure(b, 'azel') 
           # qa.angle(me.getvalue(me.measure(b,'azel'))['m0'])     # show as angles   
           # print qa.angle(me.getvalue(me.measure(b,'azel'))['m1'])  
            
          ###  az_el['m0']['value'] = (az_el['m0']['value'])+(ant_az_offsets[antenna[scan]]/360.)*2*np.pi + \
          ###    (time_az_offsets[scan]/360.)*2*np.pi
          ###    az_el['m1']['value'] = (az_el['m1']['value'])+(ant_el_offsets[antenna[scan]]/360.)*2*np.pi + \dataRA[scadataRA[scan]n]
          ###  (time_el_offsets[scan]/360.)*2*np.pi
          ###   b2 = me.measure(az_el, 'J2000')
          ###   RA_offset=dataRA[scan]-b2['m0']['value']
          ###  Dec_offset = dataDec[scan] - b2['m1']['value']
            dataRA[scan] +=(ra_offsets[scan])  #  b2['m0']['value']
            dataDec[scan] +=(dec_offsets[scan])#b2['m1']['value']
           # az_tot+=(ant_az_offsets[antenna[scan]]/360.)*2*np.pi
           # el_tot+=(ant_el_offsets[antenna[scan]]/360.)*2*np.pi
           # RA_tot += RA_offset
           # Dec_tot += Dec_offset
            ##RA_offsets = np.append(RA_offsets,RA_offset)
            ##Dec_offsets = np.append(Dec_offsets, Dec_offset)
        print ra_offsets.shape    
        plt.clf()
        plt.scatter(ra_offsets, dec_offsets, color = '#ac0c4b', linewidth = 0.0)
        plt.axhline(0, color = 'black')
        plt.axvline(0, color = 'black')
        plt.scatter(np.mean(ra_offsets),np.mean(dec_offsets) )   
        plt.ylabel('RA /degrees')
        plt.xlabel('Dec /degrees')
        plt.savefig(dirname+'RA_Dec_offsets.png') 
        
        nbins = 20
        plt.clf()
        plt.hist(ra_offsets, nbins)
        plt.savefig(dirname+'ra_offsets.png')
        plt.clf()
        plt.hist(dec_offsets, nbins)
        plt.savefig(dirname+'dec_offsets.png')
            
        if plot_azel:
            do_plot_azel(azimuths, elevations)
        # arrange back into correct array shape    
        data2 = np.dstack((np.array([dataRA]), np.array([dataDec])))
        data2 = np.transpose(data2, (2,0,1))
        print '*** data'
        print data
        print data2
        tb.putcol('DIRECTION', data2)
        tb.flush()
        
    sm.openfromms(projectname+config_file_prefix+'ms/')
    ##sm.setvp(dovp=True, usedefaultvp=True)
    # airy disk vp 
    
    #sm.setvp( dovp = True, usedefaultvp = False, vptable = output_table, dosquint = False)
    vp.setpbairy(telescope='SKA-MID',dishdiam='15m',maxrad='5deg')
    sm.setvp(dovp=True, usedefaultvp=False)
    sm.predict(complist = dirname+'points_grid.cl', incremental=False)
    print 'predicting'
    #sm.setnoise(simplenoise='1Jy')
    # sm.setgain(interval='100s', amplitude=0.01)
    sm.corrupt()
    sm.close()
    
    
def analyse(projectdir,projectname, config_file_prefix,primary_beam_size,incell_size, dirname):    
    print 'PB to first null', int(2*primary_beam_size/(incell_size))
    # perform the data reduction
    simanalyze(project=projectname, 
                        image=True, 
                        vis=projectname+config_file_prefix+'ms',  
                       # imsize = [6000,6000],
                        imsize=[int(1.2*primary_beam_size/(incell_size)),int(1.2*primary_beam_size/(incell_size))], 
                        interactive=False, 
                        niter=0, 
                        weighting='natural',  
                        pbcor= False, 
                        cell = str(incell_size)+'arcsec',
                        stokes='I', 
                        analyze=True, 
                        showuv=False, 
                        showpsf=True, 
                        showmodel=True, 
                        showconvolved=True, 
                        showclean=True, 
                        showresidual=True, 
                        showdifference=True, 
                        showfidelity=False, 
                        graphics='both', 
                        verbose = True,
                        overwrite=True)
    print projectname+config_file_prefix+'image.flat'

    os.system('rm '+'convolved_sky_model.fits')
    exportfits(imagename = projectname+config_file_prefix+'image.flat', fitsimage=projectname+'analysed_image.fits')           
    exportfits(imagename = projectname+config_file_prefix+'compskymodel.flat.regrid.conv', fitsimage='convolved_sky_model.fits')   
    
def simulator_analyse(projectname, config_file_prefix, primary_beam_size, incell_size, NUM_SPW, fieldDir, dirname, layout):
    output_ms = projectname+config_file_prefix+'ms'
    output_im = projectname+config_file_prefix
    spwmax = NUM_SPW
    # prepare input measurement set.
    im.open( output_ms )
    fieldsize = 60 #arcsec
    NPIX = int(0.25*primary_beam_size/(incell_size))
    #NPIX = int(fieldsize/incell_size)
    print NPIX
    print incell_size
    print fieldDir

    fieldDirm0 = np.copy(fieldDir['m0']['value'])
    fieldDirm1 = np.copy(fieldDir['m1']['value'])
    if layout == 'single':
        tb.open(dirname+'points_grid.cl')
        moddir= tb.getcol('Reference_Direction')
        fieldDir['m0']['value'] = np.float(moddir[0][0])
        fieldDir['m1']['value'] = np.float(moddir[1][0])
        NPIX = 512
        print fieldDir
    print fieldDir
    im.defineimage( nx = NPIX, ny = NPIX, cellx = str(incell_size)+'arcsec', celly = str(incell_size)+'arcsec', stokes = 'I', spw =   range(0,spwmax), phasecenter = fieldDir, mode = 'mfs' )
    fieldDir['m0']['value'] = np.float(fieldDirm0)
    fieldDir['m1']['value'] = np.float(fieldDirm1)
    print fieldDir
    im.setoptions( padding = 2.0 )
    im.selectvis( spw = range(0,spwmax) )
    DO_NAT = 0
    DO_UNIF= 1
    DO_ROB = 0
    if (DO_NAT):
      im.weight( type = 'natural')
    FOV = '100arcsec'
    BMAJ =  '0.8arcsec'  
    if (DO_UNIF):
      im.weight( type = 'uniform')
     # im.filter( type = 'gaussian', bmaj = BMAJ, bmin = BMAJ) 
    if (DO_ROB):
      im.weight( type = 'briggs', rmode = 'norm', robust = 0.)
    # overwrite existing images:
    os.system('rm -rf ' + output_im + 'map')
    os.system('rm -rf ' + output_im + 'map.fits')
    os.system('rm -rf ' + output_im + 'psf')
    os.system('rm -rf ' + output_im + 'psf.fits')
     
    im.makeimage( type = 'corrected', image=output_im + 'map' )
    #### error Did not get the position of SKA1-MID from data repository Frequency conversion will not work
    #### probably only a problem for spectral line imaging?
    #### if needed, might be able to input with following or similar:
    ####note, in  location
    ###For some unusual types of image, one needs to know the location to be used in calculating phase rotations. For example, one can specify images to be constructed in azel, in which case, an antenna position must be chosen. One can use functions of measures: 
  #  im.makeimage( type = 'psf', image=output_im + 'psf' )
   # params = im.fitpsf( output_im + 'psf' )
   # print 'beam fit params: ', params
    #im.makeimage( type = 'pb', image=output_im + 'pb' )
      
    im.done()
   ### im.open(output_ms)
   ### os.system('rm -r pb.map.fits')
   ### im.setvp( dovp =True, usedefaultvp = False, vptable = dirname+'vis_setup.ska1mid.-vp.tab/', dosquint = False )
   ### im.defineimage( nx = NPIX, ny = NPIX, cellx = str(incell_size)+'arcsec', celly = str(incell_size)+'arcsec', stokes = 'I', phasecenter = fieldDir)
   ### im.makeimage( type = 'pb', image='pb.map' ) 
   ### im.done()
    ia.open(output_im + 'map')
    stats = ia.statistics()
    print stats["rms"]
    ia.close()
  ###  exportfits(imagename ='pb.map', fitsimage='pb.map.fits')
    exportfits(imagename = output_im + 'map', fitsimage=output_im + 'map.fits')
   # exportfits(imagename = output_im + 'psf', fitsimage=output_im + 'psf.fits')    
 ###   allorigpb = fits.getdata('pb.map.fits')
  ###  print allorigpb
  ###  origpb = fits.getdata('pb.map.fits')[0][0]
  ###  plt.clf()
  ###  plt.imshow(origpb, cmap = 'gray')
  ###  plt.colorbar()
  ###  plt.savefig('pbarr.png')
  ###  print 'pb plotted'
        
        
    
    
def do_res(offset_value, NUM_SPW, fieldDir, primary_beam_size, incell_size, dirname, layout):    
    print 'offset: ', offset_value
    delete = 1
    observes = 1
    make_source_list = 1
    corrupts = 1
    PEs = 1
    static = 30
    time_var = 5
    plot_azel = 0
    remove_mean=0
    analyses = 1

    projectname1 = dirname+'ideal'
    projectname2 = dirname+'corrupted'
    
    os.system('cp -r '+dirname+'vis_setup'+config_file_prefix+'ms '+projectname2+config_file_prefix+'ms')    
    corrupt(projectname2, config_file_prefix, offset_value,dirname, PEs = True,  plot_azel = False)
    tb.open(projectname2+config_file_prefix+'ms/POINTING', nomodify = False)# just in case
    simulator_analyse(projectname2, config_file_prefix, primary_beam_size, incell_size ,NUM_SPW, fieldDir, dirname, layout)
    
    #get diff image
    idealfits = fits.getdata(dirname+'vis_setup'+config_file_prefix+'map.fits')
    corruptedfits = fits.getdata(projectname2+config_file_prefix+'map.fits')
    diff_image = idealfits - corruptedfits
    os.system('rm -f '+dirname+'diff_image.fits')
    fits.writeto(dirname+'diff_image.fits', diff_image)
    
    #get residuals and image
    vis_ideal = dirname+'vis_setup'+config_file_prefix+'ms'
    tb.open(vis_ideal, nomodify = True)
    ideal_corrected = tb.getcol('CORRECTED_DATA')
    tb.close()
    vis_corrupted =projectname2+config_file_prefix+'ms'
    tb.open(vis_corrupted, nomodify = False)
    corrupt_corrected = tb.getcol('CORRECTED_DATA')
    new_ideal_corrected = ideal_corrected-corrupt_corrected
    print '****'
    print ideal_corrected
    print corrupt_corrected
    print new_ideal_corrected
    tb.putcol('CORRECTED_DATA', new_ideal_corrected)
    tb.flush()
    tb.close()
    os.system('cp '+projectname2+config_file_prefix+'map.fits '+dirname+'corruptsave'+config_file_prefix+'map.fits')
    os.system('rm -f '+projectname2+config_file_prefix+'map.fits')
    
    simulator_analyse(projectname2, config_file_prefix, primary_beam_size, incell_size ,NUM_SPW, fieldDir, dirname, layout)
    res = fits.getdata(projectname2+config_file_prefix+'map.fits')[0][0]
    rms = np.std(res)
    medianabs = np.median(np.abs(res))
    maxabs = np.max(np.abs(res))
    maxx = (np.max(res))
    minn = np.min(res)
    print 'rms,maxabs: ', rms, maxabs, minn, maxx
    plt.clf()
    plt.imshow(res,cmap = 'gray_r')
    plt.colorbar()
    plt.savefig(dirname+'residuals_'+str(offset_value)+'_arcsec_normal.png')
    
    
   # os.system('ds9 ideal/idealanalysed_image.fits corrupted/corruptedanalysed_image.fits diff_image.fits' )
    return rms,maxabs, medianabs
    
i = 0
while os.path.exists("simulation%s" %i):
    i += 1
dirname =  'simulation%s/' %i
os.system('rm -rf '+dirname)
os.system('mkdir '+dirname)    
frequency = 1.4e9
dish_diameter = 15. # read this from conf file, assuming homogeneous
c = 3e8
rad2arcsec = 206265.
primary_beam_size =  (1.22*c)/(dish_diameter*frequency) # angle to first null from centre !get this from .conf
# 1.22 for circular
primary_beam_size*=rad2arcsec
print 'PB, deg: ', primary_beam_size/3600.
incell_size =0.0961 #  change to 0.09 again and change npix to 0.2*
NUM_SPW = 1  
layout = 'single'
ms_direction = 'J2000 8h00m00.031s -40d59m59.6s'
fieldDir = me.direction( ms_direction.split()[0], ms_direction.split()[1], ms_direction.split()[2])# prepare to convert to degrees 
projectname = dirname+'vis_setup' 
observe_simulator(projectname, config_file,NUM_SPW, fieldDir , dirname , layout) 
simulator_analyse(projectname, config_file_prefix, primary_beam_size, incell_size ,NUM_SPW, fieldDir, dirname, layout)  
offset_values =np.array([3000])#np.arange(10,100, 10 )  #small number problems for offset ~<5
errs = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0,64.0, 128.0, 256.0]
offset_values = errs
rmsarr = np.array([])
maxabsarr = ([]) 
medianabsarr = np.array([])
for i in offset_values:
    rms, maxabs, medianabs = do_res(i, NUM_SPW, fieldDir, primary_beam_size, incell_size, dirname, layout)
    rmsarr = np.append(rmsarr,rms)
    maxabsarr = np.append(maxabsarr, maxabs)
    medianabsarr = np.append(medianabsarr, medianabs)
print offset_values, rmsarr, maxabsarr    
plt.clf()

plt.loglog(offset_values, maxabsarr, label = 'maxabs')
plt.loglog(offset_values, medianabsarr, label = 'medianabs')
plt.loglog(offset_values, rmsarr, label = 'rms')
plt.xlim(1, 256)
plt.legend()
plt.xlabel('sigma offset (arcsec)')
plt.ylabel('error (Jy)')
plt.savefig(dirname+'plot.png')    
    
    
     
    
    
    
    
