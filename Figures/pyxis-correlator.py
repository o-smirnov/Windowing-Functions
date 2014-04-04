import Pyxis
import pyrap.tables
from Pyxis.ModSupport import *
import mqt,lsm,pyfits
# this has some useful templates we can borrow
import imager,ms
import MSResampler
#import Test
#import Source
import os
import sys
import scipy.special

#imager.DIRTY_IMAGE_Template = "${OUTFILE}-${METHOD}.dirty.fits"
DUMPFIG_Template = "${OUTFILE}-${METHOD}.wf.png"
  
# dish size, for msinfo() below
DISHSIZE = 25

def msinfo (msname=None):
  """Gets various MS info.
  Return values: NUM_TIMESLOTS,NUM_CHANNELS,MAX_BASELINE,TOTAL_TIME,TIMESLOT_SIZE,FIRST_TIMESTAMP 
  """;
  msname = msname or MS;
  tab = ms.ms(msname);
  chanwidths = ms.ms(msname,"SPECTRAL_WINDOW").getcol("CHAN_WIDTH",0,1);
  nchan = chanwidths.size;
  times = sorted(set(tab.getcol("TIME")));
  dt = times[1]-times[0];
  tottime = times[-1]-times[0]+dt;
  ntime = len(times);
  uvw = tab.getcol("UVW");
  maxbl = math.sqrt((uvw[:,:2]**2).sum(1).max());
  bw = chanwidths.ravel()[0];
  info("MS $msname: $ntime x %.2fs timeslots (%.2fh synthesis), $nchan x %.2f kHz channels, max baseline %.3f km"%
      (dt,tottime/3600,bw,maxbl/1000))
  # work out critical sampling
  # max baseline sweeps out a circle of length pi*d over 24h
  arclen = math.pi*maxbl*(tottime/(24.*3600));
  # critical sampling is at half the dish size
  nsamp = arclen/(DISHSIZE/2)
  # corresponding sampling interval
  critint = tottime/nsamp;
  info("   (NB: critical sampling interval for this baseline is %.2fs)"%critint)
  return ntime,nchan,maxbl,tottime,dt,times[0]-dt/2;

def resample_tf (mshi=None,nfreq=None,chan0=None,fov=2.,overlap=10,oper="sinc",column="DATA",outcolumn="CORRECTED_DATA",image=True,
                dump=None):
  """Resamples mshi (high-resolution MS) into MS (low-resolution MS), using normal averaging (oper="avg"),
  or one of the windowing functions defined below (oper="sinc","sinc2","airy", etc.)

  Output MS must have just 1 channel (for now). By default this all channels in the input MS, 
  but set nfreq=N to only use a block of N channels. Default then is to use the N middle
  channels, but chan0 can be used to override the starting channel.

  'data' is the input column (in mshi), 'outcolumn' is the output column (in MS).

  If 'image' is True, autoamtically makes dirty image of output MS after resampling. You can also set 'image'
  to an explicit image filename, otherwise a name is picked automatically.

  Set 'dump' to a filename to dump debugging output from the resampler.

  """;
  mshi = mshi or MSHI;
  ntime0,nfreq0,maxbl0,tottime0,dt0,time0 = msinfo(mshi);
  ntime1,nfreq1,maxbl1,tottime1,dt1,time1 = msinfo(MS);
  # support only one output frequency for now
  if nfreq1 != 1:
    abort("Output MS $MS must have one frequency (for now)")
  # averaging time factor: ratio of lo-res timeslot to hi-res timeslot
  dtime = int(dt1/dt0);
  # how many hi-res timeslots to actually use, this can be less than ntime0
  ntime = ntime1*dtime; 
  # this is the first hi-res timeslot that we start averaging at
  first_slot = int((time1-time0)/dt0);
  if first_slot < 0:
    abort("Starting time of $MS is before $mshi")
  print ntime0,first_slot,ntime;
  if first_slot+ntime > ntime0:
    abort("Ending time of $MS is after $mshi")
  info("resampling $ntime of $ntime0 timeslots by factor $dtime, starting at timeslot $first_slot");
  if nfreq is not None:
    # take center of band, if starting channel not explicitly specified
    if chan0 is None:
      chan0 = (nfreq0 - nfreq)/2;
  else:
    chan0 = 0;
    nfreq = nfreq0;
  dfreq = int(nfreq/nfreq1)
  info("resampling $nfreq channels by factor $dfreq, starting at channel $chan0")
  
  msres = MSResampler.MSResampler(mshi,column=column,time0=first_slot,ntime=ntime,freq0=chan0,nfreq=nfreq);

  if dump:
    info("Dumping debug output to file $dump");
    dump = file(dump,"w")

  if oper == "avg":
    v.METHOD = II("$column-avg-${nfreq}ch");
    info("using standard boxcar averaging");
    result = msres.boxcar(dtime,dfreq);
  elif oper in set(("sinc","sinc2","airy")):
    v.METHOD = II("$column-$oper-f${fov}-o${overlap}-${nfreq}ch");
    # define the various sinc variations
    def sinc (x,y,deg=(fov*np.pi**2)/(180.)):
      x1 = np.sqrt(x**2+y**2)*deg;
      wf = np.sin(x1)/x1;
      wf[x1==0] = 1;
      return wf;
    def sinc2 (x,y,deg=(fov*np.pi**2)/(180.)):
      x1,y1 = x*deg,y*deg;
      wx,wy = (np.sin(x1)/x1),(np.sin(y1)/y1);
      wx[x1==0] = 1;
      wy[y1==0] = 1;
      return wx*wy;
    def airy (x,y,a=2*np.pi*fov*np.pi/180):
      r = np.sqrt(x**2+y**2);
      w = scipy.special.j1(r*a)/(r*a);
      w[r==0] = 0.5;
      return w;

    func = locals().get(oper);
    if not func:
      abort("unknown oper=$oper");
    # run the actual averaging
    info("using $oper, overlap $overlap, FoV $fov");
    if overlap:
      result = msres.overlap_window(func,dtime,dfreq,overlap_time=int(dtime*overlap),dump=dump,dumpfig=DUMPFIG);
    else:
      result = msres.window(func,dtime,dfreq,dump=dump,dumpfig=DUMPFIG);
    # overlap_time=overlap*dtime,time0=first_slot,freq0=first_channel)
  else:
    abort("unknown oper=$oper");
  # save to low-res MS
  outcolumn = outcolumn or column;
  info("input column is $column, output is $outcolumn")
  MSResampler.save_visibility_arrays(v.MS,result,column=outcolumn);

  # make an image
  if image:
    imager.make_image(column=outcolumn) if type(image) is not str else imager.make_image(column=outcolumn,dirty_image=image);


#
# variables that control an resample_tf_series run
#
OPERS = ["airy","sinc2"];
FOVS = [3,2,1.5,1];
OVERLAPS = [0,1,2,4];
# WEIGHTS = "briggs",;
WEIGHTS = "natural","briggs","uniform";
RESTORE = True
# columns for which we'll make & deconvolve test images
IMAGECOLS = [ "NINESRC_DATA","UNITY_DATA","UNITY1_DATA" ];
# set this to NOISE_DATA to also make a noise image
NOISECOL = "NOISE_DATA"
## possible other columns (just as a note to myself, because that's what I simulated in the 1h9s MS):
# DATA: fine grid
# ALT_DATA: random LSM
# UNITY_DATA: one source at center
# NINESRC_DATA: 9 sources at 1.5 degrees, center source is .1 Jy.

def resample_tf_series (clearstats=False,**kw):
  """Runs a series of resampling operations using various combinations of parameters (see lists above:
  OPERS, FOVS, OVERLAPS, WEIGHTS, IMAGECOLS
  """

  v.STATSFILE = II("${DESTDIR}/stats.py");
  if clearstats:
    _clearstats();
  v.LOG = II("${DESTDIR}/log.txt");
  imager.LWIMAGER_PATH = "lwimager-trunk";
  nfreq = kw.get('nfreq');
    
  for oper in OPERS:
    for fov in (FOVS if oper != "avg" else [1]):
      for overlap in (OVERLAPS if oper != "avg" else [1]):
        # run the resampler on the imaging columns
        for col in IMAGECOLS:
          resample_tf(oper=oper,fov=fov,overlap=overlap,image=False,column=col,outcolumn="CORRECTED_DATA",**kw);
          for weight in WEIGHTS:
            imager.BASENAME_IMAGE = II("$OUTFILE-$col-$oper-f$fov-o$overlap-${nfreq>ch-}${weight}");
            imager.make_image(column="CORRECTED_DATA",weight=weight,dirty=True,restore=RESTORE); # ,npix=1024,cellsize="0.5arcsec");
        # run the resampler on the noise
        if NOISECOL:
          resample_tf(oper=oper,fov=fov,overlap=overlap,image=False,column=NOISECOL,outcolumn="MODEL_DATA",**kw);
          for weight in WEIGHTS:
            imager.BASENAME_IMAGE = II("$OUTFILE-NOISE-$oper-f$fov-o$overlap-${nfreq>ch-}${weight}");
            imager.make_image(column="MODEL_DATA",weight=weight);
            noise = pyfits.open(imager.DIRTY_IMAGE)[0].data.std();
            _writestat("noise",noise,oper,fov,overlap,weight);
            info(">>> $oper fov $fov overlap $overlap weight $weight ${nfreq <nfreq}: %.2f uJy noise"%(noise*1e+6))
      
# def run_avg_series (clearstats=False):
#   v.DESTDIR = "plots-VLAC-4h"
#   makedir(DESTDIR);
#   v.MSHI = "VLAC/VLAC-1s4h-1x10kHz-21cm.MS"
#   v.LOG = "log-VLAC.txt"
#   v.STATSFILE = "plots-VLAC-4h/stats.py"
#   v.STATQUALS = "VLAC",4
#   if clearstats:
#     _clearstats();
#   WEIGHTS = "natural","briggs"
#   for integration in 100,50,20,10: # 150,100,50,20,10:
#     # averaging
#     v.MS = makems("VLACKAT_ANTENNA","VLACa","VLAC",hours=4,integration=integration,starttime="2011/11/16/19:00");
#     average_t(column="DATA");
#     average_t(column="NOISE_DATA",outcolumn="CORRECTED_DATA");
#     imager.make_image(column="DATA",dirty_image="$DESTDIR/avg$integration-dirty.fits",wprojplanes=128);
#     noiseimg = II("$DESTDIR/avg$integration-noise.fits")
#     avgnoise = {}
#     for weight in WEIGHTS:
#       imager.make_image(column="CORRECTED_DATA",dirty_image=noiseimg,wprojplanes=0,weight=weight);
#       avgnoise[weight] = noise = pyfits.open(noiseimg)[0].data.std();
#       _writestat("noise",noise,"avg",integration,weight);
#       info(">>> averaging ($integration s) $weight weighting %.2f uJy noise"%(noise*1e+6))
  
#     # sinc filter
#     v.MS = makems("VLACKAT_ANTENNA","VLACs","VLAC",hours=4,integration=integration,starttime="2011/11/16/19:00");
#     for fov in 1,1.5,2:
#       for overlap in range(5,-1,-1):
#         average_t(column="DATA",oper="sinc",fov=fov,overlap=overlap);
#         average_t(column="NOISE_DATA",outcolumn="CORRECTED_DATA",oper="sinc",fov=fov,overlap=overlap);
#         imager.make_image(column="DATA",dirty_image="$DESTDIR/sinc$integration-F$fov-O$overlap-dirty.fits",wprojplanes=128);
#         noiseimg = II("$DESTDIR/sinc$integration-F$fov-O$overlap-noise.fits")
#         for weight in WEIGHTS:
#           imager.make_image(column="CORRECTED_DATA",dirty_image=noiseimg,wprojplanes=0,weight=weight);
#           noise = pyfits.open(noiseimg)[0].data.std();
#           _writestat("noise",noise,"sinc",integration,weight,fov,overlap);
#           _writestat("xnoise",noise/avgnoise[weight],"sinc",integration,weight,fov,overlap);
#           info(">>> sinc fov=$fov overlap=$overlap ($integration s) $weight weighting %.2f uJy noise (x%.2f)"%
#               (noise*1e+6,noise/avgnoise[weight]))

# def run_tf_series (clearstats=False):
#   v.DESTDIR = "plots-VLAC-8h"
#   makedir(DESTDIR);
#   v.MSHI = "VLAC/VLAC-1s8h-300x10kHz-21cm.MS"
#   v.LOG = "log-VLAC.txt"
#   v.STATSFILE = "plots-VLAC-4h/stats.py"
#   v.STATQUALS = "VLAC",4
#   if clearstats:
#     _clearstats();
#   WEIGHTS = "natural","briggs"
#   for integration in 100,50,20,10: # 150,100,50,20,10:
#     # averaging
#     v.MS = makems("VLACKAT_ANTENNA","VLACa","VLAC",hours=4,integration=integration,starttime="2011/11/16/19:00");
#     average_t(column="DATA");
#     average_t(column="NOISE_DATA",outcolumn="CORRECTED_DATA");
#     imager.make_image(column="DATA",dirty_image="$DESTDIR/avg$integration-dirty.fits",wprojplanes=128);
#     noiseimg = II("$DESTDIR/avg$integration-noise.fits")
#     avgnoise = {}
#     for weight in WEIGHTS:
#       imager.make_image(column="CORRECTED_DATA",dirty_image=noiseimg,wprojplanes=0,weight=weight);
#       avgnoise[weight] = noise = pyfits.open(noiseimg)[0].data.std();
#       _writestat("noise",noise,"avg",integration,weight);
#       info(">>> averaging ($integration s) $weight weighting %.2f uJy noise"%(noise*1e+6))
  
#     # sinc filter
#     v.MS = makems("VLACKAT_ANTENNA","VLACs","VLAC",hours=4,integration=integration,starttime="2011/11/16/19:00");
#     for fov in 1,1.5,2:
#       for overlap in range(5,-1,-1):
#         average_t(column="DATA",oper="sinc",fov=fov,overlap=overlap);
#         average_t(column="NOISE_DATA",outcolumn="CORRECTED_DATA",oper="sinc",fov=fov,overlap=overlap);
#         imager.make_image(column="DATA",dirty_image="$DESTDIR/sinc$integration-F$fov-O$overlap-dirty.fits",wprojplanes=128);
#         noiseimg = II("$DESTDIR/sinc$integration-F$fov-O$overlap-noise.fits")
#         for weight in WEIGHTS:
#           imager.make_image(column="CORRECTED_DATA",dirty_image=noiseimg,wprojplanes=0,weight=weight);
#           noise = pyfits.open(noiseimg)[0].data.std();
#           _writestat("noise",noise,"sinc",integration,weight,fov,overlap);
#           _writestat("xnoise",noise/avgnoise[weight],"sinc",integration,weight,fov,overlap);
#           info(">>> sinc fov=$fov overlap=$overlap ($integration s) $weight weighting %.2f uJy noise (x%.2f)"%
#               (noise*1e+6,noise/avgnoise[weight]))


define("MAKEMS_REDO",False,"if False, makems will not overwrite existing MSs");
define("MAKEMS_CONF","VLACKAT_ANTENNA","antenna table");
define("MAKEMS_NAME","VLAC","filename prefix");
define("MAKEMS_DIR","VLAC","directory for MS");

def makems (conf="$MAKEMS_CONF",name="$MAKEMS_NAME",destdir="$MAKEMS_DIR",
  hours=8,integration=60,dec=-45,freq0=1400,nchan=1,chanwidth=10,
  starttime=None):
  """Makes an MS using the specified parameters.
  hours: total synthesis duration
  integration: timeslot in seconds
  dec: declination, in degrees
  freq0: starting freq in MHz
  nchan: number of channels
  chanwidth: channel width, in kHz
  starttime: overrides default starting time (otherwise picked as transit-hours/2)
  """ 
  conf,name,destdir = interpolate_locals("conf name destdir");
  for anttab in conf,"${conf}_ANTENNA","Layouts/${conf}_ANTENNA":
    if exists(anttab):
      anttab = II(anttab);
      break;
  else:
    abort("configuration $conf not found");
  if not name:
    name = os.path.basename(anttab);
    if name.endswith("_ANTENNAS"):
      name = name.rsplit("_",1)[0];
  msname = "%s/%s-%ds%dh-%dx%dkHz-%dMHz.MS"%(destdir,name,integration,hours,nchan,chanwidth,freq0);
  info("ms $msname, configuration $conf, antenna table $anttab");
  makedir("$destdir")
  if exists(msname):
    if not MAKEMS_REDO:
      info("$msname already exists and MAKEMS_REDO=False, skipping");
      return msname;
    x.sh("rm -fr $msname");
  conffile = II("makems.${msname:BASE}.cfg");
  if not starttime:
    # work out start time: 19:00 is transit, so subtract half
    m0 = (19*60)-(hours*60)//2;
    h0 = m0//60;
    m0 = m0%60;
    starttime = "2011/11/16/%d:%02d:00"%(h0,m0);
  file(conffile,"w").write(II("""
WriteAutoCorr=Tlo
StartFreq=%g
StepFreq=%g
NFrequencies=$nchan
WriteImagingColumns=T
StepTime=$integration
#TileSizeRest=10
NParts=1
MSDesPath=.
AntennaTableName=$anttab
Declination=$dec.0.0
NBands=1
RightAscension=0:0:0
StartTime=$starttime
MSName=$msname
NTimes=%d
#TileSizeFreq=16
"""%(freq0*1e+6,chanwidth*1e+3,(hours*3600)//integration)));
  info("""creating $msname: ${hours}h synthesis, ${integration}s integration, Dec=$dec, 
  $nchan channels of $chanwidth kHz starting at $freq0 MHz""");
  # run makems
  x.makems(conffile);
  if exists(msname+"_p0") and not exists(msname):
    x.mv("${msname}_p0 $msname");
  v.MS = msname;
  # plot uv-coverage
  makedir("$DESTDIR")
  ms.plot_uvcov(ms=.1,width=10,height=10,dpi=150,save="${msname:BASE}-uvcov.png")
  return msname

def make_random_lsm (lsm="lsm.txt",nsrc=100,flux=1,dec=-45,ra=0,fieldsize=2):
  """Makes a random LSM composed of nsrc sources""";
  import numpy.random
  r = fieldsize/2
  cosdec = math.cos(dec*math.pi/180)
  xx = numpy.random.uniform(ra-r/cosdec,ra+r/cosdec,size=nsrc)
  yy = numpy.random.uniform(dec-r,dec+r,size=nsrc)

  ff = file(lsm,"w");

  ff.write("#format: ra_d dec_d i\n");
  for x,y in zip(xx,yy):
    ff.write("%f %f %f\n"%(x,y,flux));
  ff.close();
  
  info("Wrote random LSM of $nsrc sources to $lsm")



# Some useful MSs to make:
# The first one is 9 hours sampled at 1s, with 50 channels. This is the hi-res one.
#
# The others are for 8 hours (thus the high-res MS has 30 minutes overlap on each side, which allows
# overlap filters to be computed without too much bother)


# pyxis makems[hours=9,integration=1,nchan=50,chanwidth=1000]
# pyxis makems[name=VLACs,hours=8,integration=40,nchan=1,freq0=1424.5,chanwidth=1000,starttime=2011/11/16/15:00]
# pyxis makems[name=VLACa,hours=8,integration=5,nchan=1,chanwidth=1000,freq0=1424.5,starttime=2011/11/16/15:00]
# pyxis makems[name=VLACs,hours=8,integration=90,nchan=1,freq0=1424.5,chanwidth=1000,starttime=2011/11/16/15:00]
