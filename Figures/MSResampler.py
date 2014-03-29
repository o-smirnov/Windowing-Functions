import numpy as np
from pyrap.tables import table
import os
import sys
import pylab

def save_visibility_arrays (msname,arrays,column="DATA"):
  """Saves a set of visibility arrays to the MS.
  arrays is a list of (p,q,data) tuples, where p & q are antenna indices,
  and data is the corresponding array of shape ntime,nfreq,ncorr.
  The shape must match the shape of the MS, or errors will result. 
  """
  tab = table(msname,readonly=False,ack=False);
  # read in data column and antenna indices, fill data with 0
  data = tab.getcol(column)
  a1 = tab.getcol("ANTENNA1")
  a2 = tab.getcol("ANTENNA2")
  data.fill(0.);
  # fill data from arrays
  for p,q,arr in arrays:
    data[(a1==p)&(a2==q),...] = arr;
  # write out
  tab.putcol(column,data)
  tab.close();


class MSResampler (object):
    """Class for reading and resampling data from an MS"""

    def __init__ (self,msname,column="DATA",time0=0,ntime=None,freq0=0,nfreq=None):
      """Inits with given MS and column name.
      If time0/ntime and freq0/nfreq is given, handles only a subset of the data,
      otherwise reads all channels and timeslots.
      """;
      self.msname = msname;
      tab = table(msname,ack=False,readonly=False)
      self.A0 = A0 = tab.getcol('ANTENNA1')
      self.A1 = A1 = tab.getcol('ANTENNA2')
      data = tab.getcol(column);
      if nfreq is None:
        nfreq = data.shape[1]-freq0;
      self.data = data = data[:,freq0:freq0+nfreq,:];
      self.nfreq = nfreq;
      print "Visibility column shape:",data.shape
      self.na = na = np.max(A1)+1
      self.nbl = (na*(na-1))/2+na
      self.ncorr = data.shape[2]
      self.nbins = data.shape[0]
      self.UVW = tab.getcol("UVW")
      # get number of timeslots. This assumes a regular MS staructure (i.e. first baseline is
      # present)
      ntimes = (data[(A0==0)&(A1==1)]).shape[0]
      # actual number of timeslot to use given by time0 and ntime
      self.ntime = ntime or ntimes - time0;
      self.time0 = time0;

      # get frequency and wavelength (per channel) from SPECTRAL_WINDOW subtable
      t2 = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False)
      
      self.freqs = t2.getcol("CHAN_FREQ",0)[0,freq0:freq0+nfreq];
      # print "Frequencies are",self.freqs;
      self.wavel = 3e8/self.freqs;
      # print "Wavelengths are",self.wavel;
      t2.close()
      tab.close()
      
      self.F=None
      self.T='OK'

    def boxcar (self,dtime,dfreq):
      """Downsamples data using normal boxcar averaging over windows of size dtime x dfreq.
      Returns list of (antenna1,antenna2,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtime,nfreq/dfreq,ncorr)."""
      # this is the output shape
      ntime1 = self.ntime/dtime;
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # make resulting ANTENNA1/ANTENNA2 indices
      a1 = np.zeros(ntime1*self.nbl)
      a2 = np.zeros(ntime1*self.nbl)
      # loop over each baseline
      for p in range(self.na):
        for q in range(p+1,self.na):
          # extract full array for this baseline, apply subset in time
          input_index = (self.A0==p)&(self.A1==q)
          data = self.data[input_index].copy();
          data = data[self.time0:self.time0+self.ntime,...];
          # reshape so that e.g. ntime becomes ntime1,dtime, so that we can then reduce over the second axis
          data = data.reshape((ntime1,dtime,nfreq1,dfreq,self.ncorr))
          # take mean over the new axes, this becomes our result
          result.append( (p,q,data.mean(3).mean(1)) );
      return result;

    def window (self,window_function,dtime,dfreq,dump=None,dumpfig=None):
      """Downsamples data using the specified window_function. Window size is dtime x dfreq.

      window_function is a function taking an array of uv-distances (in wavelengths), and returning the
      corresponding weights.

      Returns list of (antenna1,antenna2,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtime,nfreq/dfreq,ncorr)."""
      # this is the output shape
      ntime1 = self.ntime/dtime;
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # make resulting ANTENNA1/ANTENNA2 indices
      a1 = np.zeros(ntime1*self.nbl)
      a2 = np.zeros(ntime1*self.nbl)
      # make wl: array of wavelengths over each bin -- shape x,x,NF1,DF  (where original shape was NF=NF1*DF)
      # note that the x,x stands for a "virtual" time axis
#      wl = self.wavel.reshape((nfreq1,dfreq))[np.newaxis,np.newaxis,:,:];
      wl = self.wavel.reshape((nfreq1,dfreq));
      # wl0: wavelength of centre channel of each frequency bin, shape NF1
      if dfreq%2:
        wl0 = self.wavel[dfreq/2::dfreq];
      else:
        wl0 = (self.wavel[dfreq/2-1::dfreq] + self.wavel[dfreq/2::dfreq])/2;
      # loop over each baseline
      for p in range(self.na):
        for q in range(p+1,self.na):
          # extract full array for this baseline, apply subset in time
          input_index = (self.A0==p)&(self.A1==q)
          data = self.data[input_index][self.time0:self.time0+self.ntime];
          # reshape so that e.g. NT becomes NT1,DT, so that we can then reduce over the second axis
          # resulting shape is NT1,DT,NF1,DF,NCORR
          data = data.reshape((ntime1,dtime,nfreq1,dfreq,self.ncorr))
          # get UV coordinates (omit W!) over the entire (hi-res) time range
          uv = self.UVW[input_index].copy();
          uv = uv[self.time0:self.time0+self.ntime,:2];
          # get corresponding uv-center points (uv0) for the lo-res time range, at each time window 
          # uv0 will have shape NT1,2
          if dtime%2:  # odd dtime: we have a uv-bin at the centre, so start at n/2, then take every n-th uv point
            uv0 = uv[dtime/2::dtime,:];
          else: # even dtime -- take the middle of the two centre bins
            uv0 = (uv[dtime/2-1::dtime,:] + uv[dtime/2::dtime,:])/2;
          # uv-distance has two components: along the track (in time), and across (in freq)
          # compute uvd1, the first component
          # First convert to uv-distance from window centre, along time axis, this will be of shape NT1,DT
          uvd1 = np.sqrt(((uv.reshape((ntime1,dtime,2)) - uv0[:,np.newaxis,:])**2).sum(2))
          # and now convert to wavelength at centre freq: this will be of shape NT1,DT,NF1
          uvd1 = uvd1[:,:,np.newaxis] / wl0[np.newaxis,np.newaxis,:]
          # now for the second component: ||uv0||/wl0 at each window is the centre channel uv-length,
          # the second component's distance is ||uv0||/wl - ||uv0||/wl0
          # first compute ||uv0||, shape is NT1
          uv0r = np.sqrt((uv0**2).sum(1));  
          # work out second component of uv-distance, shape will be NT1,NF1,DF
          uvd2 = uv0r[:,np.newaxis,np.newaxis]/wl[np.newaxis,:,:] - uv0r[:,np.newaxis,np.newaxis]/wl0[np.newaxis,:,np.newaxis];
          # evaluate windowing function
          # shapes are: uvd1 is NT1,DT,DF and uvd2 is NT1,NF1,DF, so insert new axes appropriately
          # wf shape is then NT1,DT,NF1,DF
          wf = window_function(uvd1[:,:,:,np.newaxis],uvd2[:,np.newaxis,:,:]);
          if dump:
            import cPickle
            cPickle.dump((p,q,uvd1[0],uvd2[0],wf[0]),dump);
          if dumpfig and p==17 and q==26:
            print uvd1[0],uvd2[0],wf[0];
            pylab.plot(wf[0,dtime/2,0,:])
            pylab.plot(wf[0,:,0,dfreq/2])
            pylab.savefig(dumpfig);
            print "Saved plot to",dumpfig;
          # sum of wf over each DTxDF bin
          wfsum = wf.sum(3).sum(1);
#          print p,q,wfsum[0]
          # apply wf to data (np.newaxis ensures every correlation is multiplied by the same wf),
          # then sum over each DTxDF bin, and normalize
          data = (data*wf[...,np.newaxis]).sum(3).sum(1)/wfsum[...,np.newaxis];
          # return result
          result.append((p,q,data));
      return result;

    def overlap_window (self,window_function,dtime,dfreq,overlap_time=0,overlap_freq=0,dump=None,dumpfig=None):
      """Downsamples data using the specified window_function. Window size is dtime x dfreq nominally,
      but overlap_time makes the effective window extend by that number of timeslots (hi-res) in each
      direction. Overlap_freq does the same for frequency (currently not implemented).

      window_function(x,y) is a function taking two array of uv-distances (in wavelengths), and returning the
      corresponding weights.

      Returns list of (antenna1,antenna2,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtime,nfreq/dfreq,ncorr)."""
      if overlap_freq:
        raise RuntimeError,"frequency overlap not yet implemented";
      # effective window size in time
      dtime1 = dtime + overlap_time*2;
      # this is the output shape
      ntime1 = self.ntime/dtime;
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # make resulting ANTENNA1/ANTENNA2 indices
      a1 = np.zeros(ntime1*self.nbl)
      a2 = np.zeros(ntime1*self.nbl)
      # make wl: array of wavelengths over each bin -- NF1,DF  (where original shape was NF=NF1*DF)
      wl = self.wavel.reshape((nfreq1,dfreq));
      # wl0: wavelength of centre channel of each frequency bin, shape NF1
      if dfreq%2:
        wl0 = self.wavel[dfreq/2::dfreq];
      else:
        wl0 = (self.wavel[dfreq/2-1::dfreq] + self.wavel[dfreq/2::dfreq])/2;
      # loop over each baseline
      for p in range(self.na):
        for q in range(p+1,self.na):
          input_index = (self.A0==p)&(self.A1==q)
          data1 = self.data[input_index].copy();
          uv1 = self.UVW[input_index,:2].copy();
          # prepare output array
	  reswind = np.zeros((ntime1,150),dtype=float)
          resdata = np.zeros((ntime1,nfreq1,self.ncorr),self.data.dtype);
          # loop over all output timeslots
          for out_time in range(ntime1):
            # make timeslice corresponding to extended window
            in_time0 = self.time0+out_time*dtime-overlap_time;
            timeslice = slice(in_time0,in_time0+dtime1);
            # extract data subset over the window
            data = data1[timeslice];
            # reshape so that NF becomes NF1,DF, so that we can then reduce over the second axis
            # resulting shape is DT1,NF1,DF,NCORR
#            print timeslice,data.shape,(dtime1,nfreq1,dfreq,self.ncorr)
            data = data.reshape((dtime1,nfreq1,dfreq,self.ncorr))
            # get UV coordinates for the window, shape is DT1,2
            uv = uv1[timeslice];
            # get corresponding uv-center point (uv0), shape is 2  
            if dtime1%2:  # odd dtime: we have a uv-bin at the centre, so start at n/2, then take every n-th uv point
              uv0 = uv[dtime1/2];
            else: # even dtime -- take the middle of the two centre bins
              uv0 = (uv[dtime1/2-1,:] + uv[dtime1/2,:])/2;
            # take uv-distance component along time, shape is DT1
            uvd1 = np.sqrt(((uv-uv0[np.newaxis,:])**2).sum(1));
            # convert it to wavelengths at centre channel of each window, this will now be of shape DT1,NF1
            uvd1 = uvd1[:,np.newaxis] / wl0[np.newaxis,:];
            # now for the second component: ||uv0||/wl0 at each window is the centre channel uv-length,
            # the second component's distance is ||uv0||/wl - ||uv0||/wl0
            # first compute ||uv0||, this is scalar
            uv0r = np.sqrt((uv0**2).sum());  
            # work out second component of uv-distance, shape will be NF1,DF
            uvd2 = uv0r/wl - uv0r/wl0[:,np.newaxis];
            # evaluate windowing function
            # shapes are: uvd1 is DT1,NF1 and uvd2 is NF1,DF so insert new axes appropriately
            # wf shape is then DT1,NF1,DF
            wf = window_function(uvd1[:,:,np.newaxis],uvd2[np.newaxis,:,:]);
            if dump and out_time==0:
              import cPickle
              cPickle.dump((p,q,uvd1,uvd2,wf),dump);
            # sum of wf over whole window and each frequency bin, shape NF1
            wfsum = wf.sum(2).sum(0);
#            if out_time == 0:
#              print p,q,wfsum;
            if dumpfig and p==17 and q==26 and out_time==0:
#             fig = pylab.figure()
#             ax = fig.add_subplot(111, projection='3d')
#             ax.plot_wireframe(wf)
              pylab.plot(wf[dtime1/2,0,:]);
              pylab.plot(wf[:,0,dfreq/2]);
              pylab.savefig(dumpfig);
              print "Saved figure",dumpfig;
              print uvd1,uvd2,wf;
            # apply wf to data (np.newaxis ensures every correlation is multiplied by the same wf),
            # then sum over whole window and each frequency bin, normalize, and assign to appropriate
            # timeslot in output
	    print wf[:,0,dfreq/2].shape
	    uvd1=uvd1[:,:,np.newaxis].reshape(uvd1.shape[0])
	    print uvd1.shape
	    if p==17 and q==26:
	#	pylab.plot(wf[dtime1/2,0,:])
		pylab.plot(uvd1,wf[:,0,dfreq/2])
		pylab.show()
	    	reswind[out_time]=wf[:,0,dfreq/2]
		result.append((p,q,wf[:,0,dfreq/2]))
            resdata[out_time,...] = (data*wf[...,np.newaxis]).sum(2).sum(0)/wfsum[...,np.newaxis];
          # return result
          #result.append((p,q,reswind));
      return result;

