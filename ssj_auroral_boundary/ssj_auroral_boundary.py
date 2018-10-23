# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # Use the non-GUI backend
from matplotlib import pyplot as pp
import matplotlib.transforms as mtransforms
from spacepy import pycdf
import logging,sys,subprocess,datetime,os,argparse,shutil,traceback
from geospacepy import special_datetime, satplottools, dmsp_spectrogram
try:
	import seaborn as sns
except:
	pass
from colorlog import ColoredFormatter #This is from PyPI, can get via `pip install colorlog`
	
#Code modified from example in Python Documentation's Logging Cookbook: 
loggername = 'ssj_auroral_boundary'
log = logging.getLogger(loggername)
log.setLevel(logging.DEBUG)

class absatday_csv(object):
	"""Class for writing satellite day of boundary identifications to CSV file
	"""
	def __init__(self,csvdir,csvfn,ssjcdffn,
					csvvars=['mlat','mlt'],
					writecsv=True):
		"""Write output CSV file

		Parameters
		----------
		csvdir : str
			Output directory (must exist)
		csvfn : str
			Output CSV filename (must exist)
		ssjcdffn : str
			DMSP SSJ CDF filename we are create boundaries for (used in header)
		csvvars : list
			List of variables to write on each line, for each boundary
			(i.e. if mlat, 4 columns will be added per row:
			(mlat of equatorward crossing 1, mlat of polar crossing 1, 
			mlat of polar crossing 2, mlat of equatorward crossing 2)
		writecsv : bool
			Calls to add_auroral_boundary_to_csv actually modify a file or not
			Default (True)
		"""
		self.csvfn = os.path.join(csvdir,csvfn)
		self.csvvars = csvvars
		self.ssjcdffn = ssjcdffn
		self.writecsv = writecsv
		if self.writecsv:
			self.write_header()

	def write_header(self):
		#Write CSV file header
		
		with open(self.csvfn,'w') as f:
			#Description
			cdffn = os.path.split(self.ssjcdffn)[-1]
			header_lines = [
				'DMSP SSJ Auroral Boundary Identification (%s)'%(cdffn),
				'Generated on %s' % (datetime.datetime.now().strftime('%c')),
				'Glossary:',
				'EQ1: First equator-side auroral boundary ',
				'PO1: First pole-side auroral boundary ',
				'PO2: Second pole-side auroral boundary ',
				'EQ2: Second equator-side auroral boundary ',
				'FOM: Figure of merit / Quality ( < 1.8 usually questionable)',
				'ind: row of DMSP SSJ CDF file corresponding to boundary',
				'hemisphere : -1 for southern hemisphere +1 for northern']
			for line in header_lines:
				f.write('# '+line+'\n')

			#Always added columns
			colnames = ('UTSecond Pass Start,'
						+'UTSecond Pass End,'
						+'hemisphere,'
						+'UTSec EQ1,'
						+'UTSec PO1,'
						+'UTSec PO2,'
						+'UTSec EQ2,'
						+'FOM,')
			#Optional columns defined in self.csvvars
			for var in self.csvvars:
				for bnd in ['EQ1','PO1','PO2','EQ2']:
					colnames+="%s %s," % (var,bnd) 
			colnames = colnames[:-1] # remove leftover comma
			colnames += '\n' # add newline
			f.write(colnames)

	def add_auroral_boundary_to_csv(self,abpp):
		"""Add a boundary idenfication to CSV

		Parameters
		----------
		abpp : abpolarpass
			Polar pass object (successful boundary identification)
		
		"""
		if self.writecsv:
			hemicode = 1 if abpp.hemi=='N' else -1
			with open(self.csvfn,'a') as f:
				line = ("%d,%d," % (abpp['uts'][0],abpp['uts'][-1])
						+"%d," % (hemicode)
						+"%d," % (int(abpp['uts'][abpp.idx_equator1]))
						+"%d," % (int(abpp['uts'][abpp.idx_pole1]))
						+"%d," % (int(abpp['uts'][abpp.idx_pole2])) 
						+"%d," % (int(abpp['uts'][abpp.idx_equator2]))
						+"%.3f" % (abpp.max_fom))
				for var in self.csvvars:
					line += ",%.3f,%.3f,%.3f,%.3f" % (
													abpp[var][abpp.idx_equator1],
													abpp[var][abpp.idx_pole1],
													abpp[var][abpp.idx_pole2],
													abpp[var][abpp.idx_equator2]
													)
				line += "\n"
				f.write(line)

class absatday(object):
	"""Class for one satellite-day of SSJ data (one CDF file)
	"""

	def __init__(self,cdffile,
					imgdir=None,make_plot=True,plot_failed=False,
					csvdir=None,writecsv=True):
		"""Class for one satellite-day of SSJ data (one CDF file)
		
		Parameters
		----------
		cdffile : str
			DMSP SSJ CDF file (probably from NASA CDAWeb)
		imgdir : str, optional
			Path to dump boundary identification images to (must exist)
			If None looks for environment variable DMSP_DIR_ABIMG
		make_plot : bool, optional
			Plot of each successful identification (the default is True)
		plot_failed : bool, optional
			Also plot unsuccesful passes (the default is False)
		writecsv : bool, optional
			Write a CSV file of boundary identifications (the default is True)
		csvdir : str, optional
			Directory to dump CSV files to (must exist)
			If None looks for environment variable DMSP_DIR_ABCSV
		"""

		self.log = logging.getLogger(loggername+'.'+self.__class__.__name__)
		self.cdf = pycdf.CDF(cdffile)
		#Parse out spacecraft so we know how to handle J4/J5 differences
		if 'dmsp' in cdffile:
			# Get the spacecraft number from the filename
			self.satnum = int(os.path.split(cdffile)[-1].split('dmsp-f')[-1][:2])
			self.log.info("Satellite number determined to be %d" % (self.satnum))
		else:
			raise RuntimeError(('Unexpected CDF filename %s, ' % (cdffile)
			                   +'could not parse out DMSP number' ))
		self.cdffn = cdffile
		self.make_plot = make_plot # Make plots of passes T/F
		self.plot_failed = plot_failed #Plot failed identifications also T/F
		self.writecsv = writecsv # Write pass identifications to a file
		self.time = self.cdf['Epoch'][:]
		self.uts = special_datetime.datetimearr2sod(self.time)
		self.hod = self.uts/3600.
		self.diff_flux = self.cdf['ELE_DIFF_ENERGY_FLUX'][:]
		self.diff_flux_std = self.cdf['ELE_DIFF_ENERGY_FLUX_STD'][:]
		self.total_flux = self.cdf['ELE_TOTAL_ENERGY_FLUX'][:]
		#The uncertainty in the CDF is relative
		self.total_flux_std = self.cdf['ELE_TOTAL_ENERGY_FLUX_STD'][:]

		#Handle filtering out any data without enough counts
		countthresh = 2.
		self.counts = (self.cdf['ELE_COUNTS_OBS'][:] \
						-self.cdf['ELE_COUNTS_BKG'][:])

		#Zero out any dubious fluxes
		self.diff_flux[self.counts <= countthresh]=0. 
		#print self.diff_flux.shape

		latvar,ltvar = 'SC_APEX_LAT','SC_APEX_MLT'
		if latvar not in self.cdf or ltvar not in self.cdf:
			#v1.1.3
			self.log.warn(('Unable to find APEX latitude or local time'
							+' variables in CDF file. Falling back to'
							+' AACGM magnetic coordinates'))
			latvar,ltvar = 'SC_AACGM_LAT','SC_AACGM_LTIME'
			 
		self.mlat = self.cdf[latvar][:]
		self.mlt = self.cdf[ltvar][:]
		self.channel_energies = self.cdf['CHANNEL_ENERGIES'][:]
		self.xings = self.simple_passes(self.mlat)
		self.polarpasses = []
		
		#Look for environemnt variables to define paths if no paths provided
		imgdir = self.if_none_use_envvar(imgdir,'DMSP_DIR_ABIMG')
		if imgdir is None:
			raise RuntimeError('No image dir passed & no DMSP_DIR_ABIMG envvar')
		self.imgdir = imgdir

		csvdir = self.if_none_use_envvar(csvdir,'DMSP_DIR_ABCSV')
		if csvdir is None:
			raise RuntimeError('No csv dir passed & no DMSP_DIR_ABCSV envvar')
		
		cdffn_noext = os.path.splitext(os.path.split(cdffile)[-1])[0]
		csvfile = cdffn_noext+'_boundaries.csv'
		self.csv = absatday_csv(csvdir,csvfile,cdffile,writecsv=self.writecsv)

		#Start processing the polar passes one by one
		for i in range(len(self.xings)-1):
			newpass = abpolarpass(self,self.xings[i],self.xings[i+1]-1)
			self.polarpasses.append(newpass)

	def if_none_use_envvar(self,checkvar,envvar):
		"""Check for environment variable envvar if checkvar is None
		"""
		if checkvar is None and envvar in os.environ:
				return os.environ[envvar]
		elif checkvar is None:
			return None
		else:
			return checkvar

	def simple_passes(self,latitude):
		"""Finds all the equator crossings"""
		npts = len(latitude.flatten())
		entered_north = []
		entered_south = []
		
		for k in range(1,npts):
			#poleward crossing
			if latitude[k-1] < 0. and latitude[k] >= 0.:
				entered_north.append(k)
				self.log.info( "Entered Northern Hemisphere: ind:%d,lat:%.3f" % (k,latitude[k]) )
			elif latitude[k-1] > 0. and latitude[k] <= 0.:
				entered_south.append(k)
				self.log.info("Entered Southern Hemisphere: ind:%d,lat:%.3f" % (k,latitude[k]))

		xings = entered_north+entered_south
		xings.sort()
		return xings

	def __getitem__(self,var):
		if hasattr(self,var):
			return getattr(self,var)
		elif var in self.cdf:
			return self.cdf[var][:]
		else:
			self.log.error(("Non-existent variable %s" % (str(var))
							+" requested through getattr. Returning None"))
			return None
		   
class absegment(object):
	"""
	A class for encapulating one 'segment', 
	that is one region of over-threshold flux
	"""
	def __init__(self,polarpass,ind_seg_start,ind_seg_end):
		self.log = logging.getLogger(loggername+'.'+self.__class__.__name__)
		self.polarpass = polarpass
		self.si = ind_seg_start
		self.ei = ind_seg_end
		self.area = np.nansum(self['intflux']) # Total flux in this segment 
		self.area_uncert = np.nanmean(self['total_flux_std']) # Avg Relative Unc
		self.twidth = abs(self['uts'][-1]-self['uts'][0]) #Width in seconds
		self.log.debug(('Segment from UT second '
						+'%.1f - %.1f ' % (self['uts'][0],self['uts'][-1])
						+'(%.1f seconds) ' % (self.twidth)
						+'had area (total flux) %.2g' % (self.area)))
		
	def __str__(self):
		return "UTS: %.1f - %.1f (DT=%.1f s), area (totflux): %.2g" % (self['uts'][0],self['uts'][-1],self.twidth,self.area)

	def __getitem__(self,var):
		"""Call the parent's getitem for the variable, and then subscript it"""
		return self.polarpass[var][self.si:self.ei]

class abpolarpass(object):
	"""A class for encapulating each spacecraft pass over the pole. 
		Handles finding the boundaries for a single pass
	"""
	def __init__(self,satday,ind_pass_start,ind_pass_end,settings=None):
		
		self.log = logging.getLogger(loggername+'.'+self.__class__.__name__)
		self.si = ind_pass_start # Index into satday's data
		self.ei = ind_pass_end #Index into satday's data
		self.satday = satday # Parent spacecraft day of data
		
		self.hemi = 'N' if np.nanmean(self['mlat'])>1 else 'S'
		self.idx_pole_approach = np.argmax(np.abs(self['mlat']))
		
		if settings is None:
			settings = dict()
		
		#Boundary Finding Settings
		self.MAX_DATA_GAP = 60 #Max number of missing seconds in top 9 channels
		#self.FLUX_MIN = 10**6.5 # Original (possibly appropriate for J4?)
		self.FLUX_MIN = 10**9  
		self.MIN_IND_GAP = 5     # minimum # of samples between individial segments
		self.MIN_SAMPLES = 45    # minimum # samples within a peak, 1 sample = 1 second
		self.MIN_LAT     = 50.   # start search
		self.MIN_DT      = 120.  # minimum seconds between end of peak[n] and start of 
							#peak[n+1].  (smaller than this is a skimmer pass)
		self.MIN_DT_AREA = 30.   # minimum seconds of an area.  
							#(smaller than this should be ignored) 
		self.MIN_EQ_SPIKE_DT = 30. # used for selecting equatorward edge of the auroral zone

		#Define empty versions of everything

		#Step 1: prepare data
		self.intflux,self.total_missing_samples,self.longest_strech_missing = None,None,None

		#Step 2: Find Segments
		self.segments,self.max_seg_area = None,None
		
		#Step 3: Compare Segments and Compute Figure of Merit
		self.foms,self.combos,self.max_fom,self.max_fom_ind = None,None,None,None

		#Step 4: Compute the actual boundaries, filtering for possible short duration spikes
		self.idx_pole1,self.idx_equator1,self.idx_pole2,self.idx_equator2 = None,None,None,None

		self.failure_reason = None # This is set in the methods then read in the plot method

		#Integrate the flux and determine how much data is available in this pass
		self.intflux,self.total_missing_samples,self.longest_strech_missing = self.prepare_flux_data()
		if self.intflux is not None:
			#Find the segments	    
			self.segments,self.max_seg_area = self.find_segments() 
			#Compare the segments and compute the figure of merit for each combo
			if self.segments is not None:
				self.foms,self.combos,self.max_fom,self.max_fom_ind = self.compare_segments()
				#If we have a good set of combos (i.e. non-nan figure of merit), go and finalize the boundaries
				if self.max_fom is not None:
					self.idx_pole1,self.idx_equator1,self.idx_pole2,self.idx_equator2,self.segment1,self.segment2 = self.get_boundaries()
				else:
					self.log.error("No finite FOM found! Pass processing for %s is discontinued" % (str(self)))
			else:
				self.log.error("Segment finding was aborted (returned None). Pass processing for %s is discontinued" % (str(self)))
				
		else:
			self.log.error("Failed to integrate flux...too many missing values")
			

		if self.intflux is not None and self.satday.make_plot:
			if self.failure_reason is None or self.satday.plot_failed:
				f= self.plot()

				self.figfile = os.path.splitext(os.path.split(self.satday.cdffn)[-1])[0]+'_%spass_uts%.5d_uts%.5d.png' % (self.hemi,np.floor(self['uts'][0]),np.floor(self['uts'][-1]))
				self.figfile = os.path.join(self.satday.imgdir,self.figfile)
				self.log.debug("Figfile is %s" % (self.figfile) )
				f.savefig(self.figfile,dpi=300.)
				pp.close(f)
		
		boundary_indices = [self.idx_equator1,
							self.idx_pole1,
							self.idx_pole2,
							self.idx_equator2]

		if self.failure_reason is None and None not in boundary_indices: 
			self.satday.csv.add_auroral_boundary_to_csv(self)

	def draw_boundary_vlines(self,ax,timevar='uts',inplace_label=True):
		"""
		Draw vertical lines at the detected boundaries on an arbitrary axes
		timevar must be a valid variable to get through the __getitem__ interface
		"""
		mtrans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
		boundary_idxs = [self.idx_equator1,self.idx_pole1,self.idx_pole2,self.idx_equator2]
		boundary_names = ['EQ1','PO1',"PO2",'EQ2']
		boundary_colors = ['b','r','g','m']
		boundary_alignments = ['right','left','right','left']
		lw=1.5
		txtshift=5
		for idx,name,color,alignment in zip(boundary_idxs,boundary_names,boundary_colors,boundary_alignments):
			if idx is not None:
				ax.axvline(self[timevar][idx],color=color,label=name,lw=lw)
				if inplace_label:
					ax.text(self[timevar][idx],1.1,name,color=color,horizontalalignment=alignment,transform=mtrans)

	def get_boundary_reports(self):
		desc1 = "Boundary 1: %.3f-%.3f, A1: %.1e, RelUncertA1: %.1f, A1/A_max: %.3f, twidth1: %.3fs" % (self['mlat'][self.idx_equator1],self['mlat'][self.idx_pole1],
																			self.segment1.area,self.segment1.area_uncert,self.segment1.area/self.max_seg_area,
																			self.segment1.twidth)
		desc2 = "Boundary 2: %.3f-%.3f, A2: %.1e, RelUncertA2: %.1f, A2/A_max: %.3f, twidth2: %.3fs" % (self['mlat'][self.idx_equator2],self['mlat'][self.idx_pole2],
																		self.segment2.area,self.segment2.area_uncert,self.segment2.area/self.max_seg_area,
																		self.segment2.twidth)
		descpc = "Polar Cap Width in Time %.1f sec, Boundary Set Score/FOM %.1f" % (self.segment2['uts'][0]-self.segment1['uts'][-1],
																								self.max_fom)
		return desc1,desc2,descpc

	def shade_boundary_canidates(self,ax,timevar='uts',color='green',alpha=.2,fs=12):
		mtrans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
		if self.segments is not None:
			for i,seg in enumerate(self.segments):
				ax.axvspan(seg[timevar][0],seg[timevar][-1],facecolor=color,
							alpha=alpha)
				ax.text((seg[timevar][-1]+seg[timevar][0])/2.,.05,str(i),
						color='black',fontsize=fs,transform=mtrans,ha='center')

	def plot(self):
		"""
		Draw a dialplot of the data
		"""
		
		f = pp.figure(figsize=(8.5,11),dpi=300)
		a = f.add_subplot(3,1,1)
		a2 = f.add_subplot(3,1,2)
		a3 = f.add_subplot(3,1,3)
		
		iflx = self['intflux']
		iflx[iflx==0.] = .1

		ms=2

		#Plot the track
		X,Y = satplottools.latlt2cart(self['mlat'],self['mlt'],self.hemi)
		a.plot(X,Y,'k.',markersize=ms)

		mappable = satplottools.hairplot(a,self['mlat'],self['mlt'],np.log10(self['intflux']),self.hemi,
			vmin=np.log10(self.FLUX_MIN),vmax=12)
		f.colorbar(mappable,label='log10(Smoothed Integrated EEFlux)',ax=a)

		titlstr = "%s\n" % (os.path.split(self.satday.cdffn)[1])

		a.text(0,0,self.hemi)
		#def dmsp_spectrogram( times, flux, channel_energies, lat=None, lt=None, fluxunits='eV/cm^2-s-sr-eV',
		#		logy=True, datalabel=None, cblims=None, title=None, ax=None, ax_cb=None ):
		
		fluxstd = self.moving_average(self['total_flux_std'],15)
		a2.plot(self['uts'],self['intflux'],'k.',label='Smoothed Int >1KeV Flx',ms=5.)
		a2.plot(self['uts'],fluxstd,'r.',ms=3,label='Relative Uncertainty')
		#a2.axhline(np.nanmean(fluxstd)-.5*np.nanstd(fluxstd),label='Mean Uncertainty',color='orange')
		#a2.plot(self['uts'],self['intflux']-3*self['intflux']*fluxstd,'g.',ms=5,label='3*std lower bound')				
		a2.set_xlabel("UT Second of Day")
		a2.set_yscale('log')
		a2.axhline(self.FLUX_MIN,label='Threshold',color='grey')
		a2.set_title("Integrated Flux (9 Highest E Channels %.2feV-%.2feV)" % (self['channel_energies'][0],self['channel_energies'][8]))

		dmsp_spectrogram.dmsp_spectrogram(self['time'],self['diff_flux'],self['channel_energies'],lat=self['mlat'],lt=self['mlt'],
			ax=a3,cblims=[1e5,1e10])
		#a3.set_ylim([1e2,1e5])	

		if self.segments is not None:
			for seg in self.segments:
				a2.axvspan(seg['uts'][0],seg['uts'][-1],facecolor='green',alpha=.2)

		a2.grid(True)
		a3.grid(True)

		lw=1.5
		txtshift=5
		#Add boundaries
		if self.idx_equator1 is not None:
			a.plot(X[self.idx_equator1],Y[self.idx_equator1],'bo',alpha=.5)
			a2.axvline(self['uts'][self.idx_equator1],color='b',label='EQ1',lw=lw)
			a3.axvline(self['time'][self.idx_equator1],color='b',label='EQ1',lw=lw)
			a.text(X[self.idx_equator1],Y[self.idx_equator1]-txtshift,"EQ1",color='b',horizontalalignment='left')
		else:
			titlstr+='No 1st EQB,'

		if self.idx_pole1 is not None:
			a.plot(X[self.idx_pole1],Y[self.idx_pole1],'ro',alpha=.5)
			a2.axvline(self['uts'][self.idx_pole1],color='r',label='PO1',lw=lw)
			a3.axvline(self['time'][self.idx_pole1],color='r',label='PO1',lw=lw)
			a.text(X[self.idx_pole1],Y[self.idx_pole1]-txtshift,"PO1",color='r',horizontalalignment='right')
		else:
			titlstr+='No 1st PWB,'

		if self.idx_pole2 is not None:
			a.plot(X[self.idx_pole2],Y[self.idx_pole2],'go',alpha=.5)
			a2.axvline(self['uts'][self.idx_pole2],color='g',label='PO2',lw=lw)
			a3.axvline(self['time'][self.idx_pole2],color='g',label='PO2',lw=lw)
			a.text(X[self.idx_pole2],Y[self.idx_pole2]-txtshift,"PO2",color='g',horizontalalignment='left')
		else:
			titlstr+='No 2nd PWB,'

		if self.idx_equator2 is not None:
			a.plot(X[self.idx_equator2],Y[self.idx_equator2],'mo',alpha=.5)
			a2.axvline(self['uts'][self.idx_equator2],color='m',label='EQ2',lw=lw)
			a3.axvline(self['time'][self.idx_equator2],color='m',label='EQ2',lw=lw)
			a.text(X[self.idx_equator2],Y[self.idx_equator2]-txtshift,"EQ2",color='m',horizontalalignment='right')
		else:
			titlstr+='No 2nd EQB,'

		if titlstr[-1]==',':
			titlstr = titlstr[:-1]+'\n' # Remove trailing comma, add newline

		a.legend(ncol=2)
		a2.legend(loc=0)

		if self.failure_reason is not None:
			a.text(-40,-60,self.failure_reason,color='red')
		else:
			desc = "Boundary 1: %.3f-%.3f, A: %.1e, RelUncertA: %.1f, A/A_max: %.3f, twidth: %.3fs\n" % (self['mlat'][self.idx_equator1],self['mlat'][self.idx_pole1],
																			self.segment1.area,self.segment1.area_uncert,self.segment1.area/self.max_seg_area,
																			self.segment1.twidth)
			desc += "Boundary 2: %.3f-%.3f,A: %.1e, RelUncertA: %.1f, A/A_max: %.3f, twidth: %.3fs\n" % (self['mlat'][self.idx_equator2],self['mlat'][self.idx_pole2],
																			self.segment2.area,self.segment2.area_uncert,self.segment2.area/self.max_seg_area,
																			self.segment2.twidth)
			desc += "Polar Cap Width in Time %.1f sec" % (self.segment2['uts'][0]-self.segment1['uts'][-1])
			a.text(-40.,-60.,desc,color='blue')

		if self.max_fom is not None:
			titlstr += 'Identification FOM ( < 1.8 is questionable ): %.2f' % (self.max_fom)
		
		f.suptitle(titlstr)

		#f.autofmt_xdate()
		f.tight_layout()

		return f

	def __str__(self):
		return "%s hemisphere pass UTS: %.3f-%.3f, MLT: %.3f-%.3f, closest pole approach mlat:%.3f @UTS:%.3f" % \
			(self.hemi,self['uts'][0],self['uts'][-1],
				self['mlt'][0],self['mlt'][-1],
				self['mlat'][self.idx_pole_approach],self['uts'][self.idx_pole_approach])

	def __getitem__(self,var):
		#Calls the parent satday's __getitem__ method for the desired variable
		#this way we have a uniform interface, and can always get the right length
		#of data without storing it in this object (memory optimization)
		#Variable is an attribute of the polar pass
		if var in ['intflux','intflux_std']:
			return getattr(self,var)
		#Variable is an attribute of satday, but not indexable by si:ei
		elif var == 'channel_energies':
			return self.satday[var] # Not the same length (not time-varying)
		#Variable is an attribute of satday and we want only the part for this pass
		else:
			return self.satday[var][self.si:self.ei+1]

	def prepare_flux_data(self):
		#Compute indices of all flux columns' bad data, determine if there are too many gaps
		channel_gaps_unacceptable = []
		badinds,goodinds = dict(),dict()
		diff_eflux = self['diff_flux']
		channel_energies = self['channel_energies']

		total_missing_samples = 0
		longest_strech_missing = 0
		for c in range(len(diff_eflux[0,:])): #Iterate over columns
			badinds[c] = np.flatnonzero(np.logical_not(np.isfinite(diff_eflux[:,c])))
			goodinds[c] = np.flatnonzero(np.isfinite(diff_eflux[:,c]))
			self.log.info('Channel #%d (E=%.3feV) had %d nan points' % (c,channel_energies[c],len(badinds[c])))
			# ;Time gaps:  Look for them in each of the 9 highest channels since those are the channels we're using
			# ;TODO: look for 1) total missing number of samples; 2) longest stretch of missing data
			channel_gaps_unacceptable.append(len(badinds[c]) > self.MAX_DATA_GAP)
			total_missing_samples += len(badinds[c])
			#Find the largest difference in time between good values
			if len(badinds[c])>1:
				longest_strech_missing = longest_strech_missing if longest_strech_missing > np.nanmax(np.diff(badinds[c])) else np.nanmax(np.diff(badinds[c]))  

		ngaps_unacceptable = all(channel_gaps_unacceptable)
		if ngaps_unacceptable:
			self.log.error('Highest 9 SSJ channels ALL had more than %d missing seconds of data, this polar pass will be ignored' % (self.MAX_DATA_GAP))
			self.failure_reason = 'Too many missing seconds (>%d) of data in SSJ channels' % (self.MAX_DATA_GAP)
			return None,None,None

		# Integrate the flux
		#Calling '__getitem__', i.e. using self as a dictionary (self['key']) returns:
		#	All columns from the parent satday variable satday.key, for only the rows that correspond to this pass
		intflux_notsmooth = self.flux_integrate(self['diff_flux'],self['channel_energies'])
		#Integrate the uncertainty
		#intflux_std_notsmooth = self.std_integrate(self['diff_flux_std'],self['channel_energies'])
		
		self.log.debug("In prepare_flux_data, shape intflux is %s" % (str(intflux_notsmooth.shape)))
		#Fail if there aren't enough points to smooth
		if len(intflux_notsmooth) < 15:
			self.log.error('Not enough points to smooth %d < 15' % (len(intflux_notsmooth)))
			self.failure_reason = 'Not enough data to smooth integrated flux!'
			return None, None, None

		intflux = self.moving_average(intflux_notsmooth,15)
		#intflux_std = self.moving_average(intflux_std_notsmooth,15)
		return intflux,total_missing_samples,longest_strech_missing

	def std_integrate(self,diff_flux_std,energies):
		"""
		Integrate the UNCERTAINTY in differential flux using Hardy algorithm for
		channel energy width.
		Just applying the standard uncertainty propagation formula
		Assumes no covariance between channel fluxes (as done with SSJ uncertainty creation)

		sJE(je_i)^2 = SUM[ (dJE/dje_i)^2*sje_i^2 ]

		where s represents uncertainty
		JE is integrated flux
		je_i is differential flux of i-th SSJ channel
		and the SUM runs over i in 0 (30KeV) to 8(1.392 KeV) 
		"""
		c = energies.copy()

		#Turn all NaN to zeros for purposes of integration
		nantozero_std = diff_flux_std.copy()
		for ch in range(len(nantozero_std[0,:])):
			bad = np.logical_not(np.isfinite(nantozero_std[:,ch]))
			neg = nantozero_std[:,ch]<0.
			
			nantozero_std[np.logical_or(bad,neg),ch] = 0.

		#Integrate
		std_integrated  = (c[0] - c[1])**2*nantozero_std[:,0]**2 +\
							(1./2*(c[0] - c[2]))**2*nantozero_std[:,1]**2 +\
							(1./2*(c[1] - c[3]))**2*nantozero_std[:,2]**2 +\
							(1./2*(c[2] - c[4]))**2*nantozero_std[:,3]**2 +\
							(1./2*(c[3] - c[5]))**2*nantozero_std[:,4]**2 +\
							(1./2*(c[4] - c[6]))**2*nantozero_std[:,5]**2 +\
							(1./2*(c[5] - c[7]))**2*nantozero_std[:,6]**2 +\
							(1./2*(c[6] - c[8]))**2*nantozero_std[:,7]**2 +\
								 (c[7] - c[8])**2*nantozero_std[:,8]**2
		return np.sqrt(std_integrated)

	def flux_integrate(self,diff_flux,energies):
		"""Integrate the differential flux using Hardy algorithm for
		channel energy width
		From Redmon IDL version:
		; Integrates UPPER (30keV to 1.39keV) Differential energy fluxes.
		; ASSUMES these channels are the 0th - 8th columns in the input matrix.
		; Does NOT do any smoothing or other funny business!
		"""
		c = energies.copy()

		#Turn all NaN to zeros for purposes of integration
		nantozero_flux = diff_flux.copy()
		for ch in range(len(nantozero_flux[0,:])):
			bad = np.logical_not(np.isfinite(nantozero_flux[:,ch]))
			neg = nantozero_flux[:,ch]<0.
			
			nantozero_flux[np.logical_or(bad,neg),ch] = 0.

		#Integrate
		flux_integrated  = (c[0] - c[1])*nantozero_flux[:,0] +\
							1./2*(c[0] - c[2])*nantozero_flux[:,1] +\
							1./2*(c[1] - c[3])*nantozero_flux[:,2] +\
							1./2*(c[2] - c[4])*nantozero_flux[:,3] +\
							1./2*(c[3] - c[5])*nantozero_flux[:,4] +\
							1./2*(c[4] - c[6])*nantozero_flux[:,5] +\
							1./2*(c[5] - c[7])*nantozero_flux[:,6] +\
							1./2*(c[6] - c[8])*nantozero_flux[:,7] +\
								 (c[7] - c[8])*nantozero_flux[:,8]
		
		return flux_integrated

	def rolling_window(self,a, window):
		shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
		strides = a.strides + (a.strides[-1],)
		return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

	def moving_average(self,x,window_size):
		"""Creates a weighted average smoothed version of x using the weights in window"""
		return np.nanmean(self.rolling_window(np.concatenate((x[:window_size/2],x,x[-window_size/2+1:])),window_size),-1)

	def find_segments(self):
		uts = self['uts']
		lat = self['mlat']
		mlt = self['mlt']
		intflux = self['intflux']
		stdflux = self.moving_average(self['total_flux_std'],15) # Smoothed 
		#lowerbndflux = intflux - 3*stdflux*intflux #3std lower bound

		#self.log.debug('In find_segments intflux.shape is %s, lat.shape is %s ' % (str(intflux.shape),str(lat.shape)))
		idx_pass    = np.flatnonzero( np.abs(lat) > self.MIN_LAT )
		abovethresh = intflux - self.FLUX_MIN > 0.
		idx_crosses = np.flatnonzero( np.logical_and( abovethresh, np.abs(lat) > self.MIN_LAT ) )

		#self.log.debug('\nidx_crosses: %s\n' % (str(idx_crosses)))
 
		# Not enough flux > threshold => SKIP
		if len(idx_crosses) < self.MIN_SAMPLES:
			self.log.info("not enough points %d > threshold %d => SKIPPING" % (len(idx_crosses),self.MIN_SAMPLES))
			self.failure_reason = "Number of points with above threshold flux (%d found < %d minimum)" % (len(idx_crosses),self.MIN_SAMPLES)
			return None,None
	
		# Determine the end of each segment by differencing the indices and seeing when the adjacent flux 
		# over threshold are more than 1 time step apart
		# each region of continuous 'over-threshold' fluxes is henceforth called a segment
		deltaflux = np.abs( np.diff(idx_crosses) )
		idx_segment_ends  = np.flatnonzero( deltaflux > self.MIN_IND_GAP )

		idx_segment_ends = np.concatenate((idx_segment_ends, np.array([len(idx_crosses)-1]))) # end of last segment is last point above threshold
		n_segments = len( idx_segment_ends )

		#start of each segment
		idx_segment_starts = np.concatenate( (np.array([0]),idx_segment_ends[:-1] + 1) ) #Start is 1 after end

		self.log.debug('\nidx_segment_starts: %s\n' % (str(idx_segment_starts)))
		self.log.debug('\nidx_segment_ends: %s\n' % (str(idx_segment_ends)))
 
		# Fail because not enough segments found
		if n_segments <= 2:
			self.log.info("not enough above threshold intervals greater than %d points wide" % (self.MIN_IND_GAP))
			self.failure_reason = "Not enough above-threshold intervals (<2) greater than %d points wide" % (self.MIN_IND_GAP)
			return None,None

		# Remove segments that start at the first latitude of the search
		if idx_crosses[ idx_segment_starts[ 0 ] ] == idx_pass[ 0 ]:
			self.log.info("first segment starts at first latitude of search => REMOVING")
			idx_segment_starts = idx_segment_starts[ 1:n_segments ]
			idx_segment_ends   = idx_segment_ends[   1:n_segments ]
			n_segments -= 1
		
		# Remove segments that end at the last latitude of the search
		if idx_crosses[ idx_segment_ends[-1] ] == idx_pass[-1]:
			self.log.info("last segment ends at last latitude of search => REMOVING")
			idx_segment_starts = idx_segment_starts[ 0:n_segments - 1 ]
			idx_segment_ends   = idx_segment_ends[   0:n_segments - 1 ]
			n_segments -= 1
		
		# Remove segments that have the same start and end, or end = start+1 (length 1 segments)
		# This has to be an iterative process because the indices change each time we resize the idx arrays
		lenzerosegs = np.flatnonzero(idx_segment_starts==idx_segment_ends)
		while len(lenzerosegs)>0:
			ind = lenzerosegs[0]
			self.log.error("Removing segment %d because start and end are the same (%d==%d)" % (ind,idx_segment_starts[ind],idx_segment_ends[ind]))
			idx_segment_starts = np.concatenate( (idx_segment_starts[ :ind ], idx_segment_starts[ ind+1: ]) )
			idx_segment_ends = np.concatenate(  (idx_segment_ends[ :ind ], idx_segment_ends[ ind+1: ] ) )
			n_segments -= 1
			lenzerosegs = np.flatnonzero(idx_segment_starts==idx_segment_ends)
			
		# Added because of F12199711070203.J4

		# 1 Segment skimmer passes => Skip
		if ( n_segments < 1 ):
			self.log.error("only found one segment => SKIPPING")
			self.failure_reason = "Only one region of above threshold flux (unibrow)"
			return None,None
		
		# "unibrow" skimmer pass (i.e. humps too close together) => Skip
		if np.max( deltaflux ) < self.MIN_DT :
			self.log.error("skimmer pass => SKIPPING")
			self.failure_reason = "Regions of above threshold flux too close.\n(dt: %.1f < %.1f) (borderline unibrow)" % (np.max(deltaflux),self.MIN_DT) 
			return None,None
		
		#If we got this far we can turn self.segments into something other than None
		#Nope moved to a pass-back sort of paradigm
		
		segments = []
		max_seg_area = 0.
		for k in range(n_segments):	
			si,ei = idx_crosses[idx_segment_starts[k]],idx_crosses[idx_segment_ends[k]]
			self.log.debug("Now processing segment %d: idx_segment_starts = %d, idx_segment_ends = %d" % (k,si,ei))
			segments.append(absegment(self,si,ei))
			self.log.debug("\n--Added Segment #%d: (t[%d]=%.1f,lat[%d]=%.3f) - (t[%d]=%.1f,lat[%d]=%.3f), delta_t = %s seconds" % (k,
				si,uts[si],si,lat[si],ei,uts[ei],ei,lat[ei],uts[ei]-uts[si]))
			new_max_seg_area = max([abs(max_seg_area),abs(segments[-1].area)])
			if max_seg_area < new_max_seg_area:
				self.log.debug("Segment %d: New maximum total flux (area, used in FOM determination) : %.3g" % (k,max_seg_area))
				max_seg_area = new_max_seg_area
			
		return segments,max_seg_area

	def compare_segments(self):
		"""Compare every possible pair of segments to determine which best pair best represents the boundaries"""
		combos = []
		foms = []
		
		#Edge case: only one segment identified, so can't compare
		if len(self.segments)==1:
			fom_failure_reason = 'Only one segment found, so could not compare'
			
		for s1_ind,s1 in enumerate(self.segments):
			for k,s2 in enumerate(self.segments[s1_ind+1:]):
				s2_ind = s1_ind+1+k #Actual index in segments of s2
				#Time from poleward edge (end) of ascending to poleward edge (beginning) of descending
				inside_twidth = s2['uts'][0]-s1['uts'][-1] 
				#Require that both segements are wider than a tuning parameter
				#This prevents isolated spikes that maximize the fom because they are low latitude
				#From being picked
				if s1.twidth < self.MIN_DT_AREA or s2.twidth < self.MIN_DT_AREA: 
					fom = np.nan
					fom_failure_reason = "Areas are too short < %s " % (str(self.MIN_DT_AREA))
				#Require that the two segments to be compared have the maximum latitude point between them
				#this prevents odd off center identifications
				elif self.idx_pole_approach <= s1.ei or self.idx_pole_approach >= s2.si:
					fom = np.nan
					fom_failure_reason = "Off center, maximum latitude location not between segements"
				else:
					# If the segments aren't too narrow in time or too close together
					fom = (s1.area + s2.area)/self.max_seg_area + \
							1.-s1.area_uncert + 1.-s2.area_uncert + \
							inside_twidth/(20.*60.) #Figure of merit for this combo
					foms.append(fom)
					combos.append((s1_ind,s2_ind))

				self.log.debug("Comparing \n#%d:%s and \n#%d:%s, \nFOM is %.3f" % (s1_ind,str(s1),s2_ind,str(s2),fom))
				if not np.isfinite(fom):
					self.log.warn(fom_failure_reason)

		foms = np.array(foms)
		if len(foms)>0:
			max_fom = np.nanmax(foms) 
			max_fom_ind = np.flatnonzero(foms==max_fom)[0]
			self.log.info("Max fom was %.3f, matching combo of segments %d and %d" % (max_fom,combos[max_fom_ind][0],combos[max_fom_ind][1]))
		else:
			max_fom = None
			max_fom_ind = None
			self.failure_reason = 'No Valid Comparisons (b/c %s)' % (fom_failure_reason)
			self.log.warn("No non-nan FOMs found! Means %s" % (fom_failure_reason))
		
		return foms,combos,max_fom,max_fom_ind

	def get_boundaries(self):
		"""The final step in getting the boundaries. Returns the indices into the pass' data for the 
		poleward and equatorward edges of the boundaries"""

		pole1segnum = self.combos[self.max_fom_ind][0]
		pole2segnum = self.combos[self.max_fom_ind][1]

		segpole1 = self.segments[pole1segnum]
		segpole2 = self.segments[pole2segnum]
		
		#Area1 => Ascending phase
		# pick eq1 ignoring very narrow spikes
		
		#Start by assuming that the equatorward boundary is
		#simply the start of the segment for which we found the first poleward boundary
		#If we don't find anything better equatorward of that we will stick with it
		idx_equator1=segpole1.si
		for i,seg in enumerate(self.segments[:pole1segnum+1]):
			self.log.info("Found first equatorward boundary at segment #%d:%s" % (i,str(seg)))
			if ( seg.twidth >= self.MIN_EQ_SPIKE_DT and seg.area_uncert < np.nanmean(self['total_flux_std'])):
				idx_equator1 = seg.si     #if meets time requirement, define as boundary
				break

		if idx_equator1 is None:
			self.log.warn('No 1st equatorward boundary found satisfying delta-t > %s' % (str(self.MIN_EQ_SPIKE_DT)))
			
		
		idx_pole1   = segpole1.ei
		self.log.info("First poleward boundary at segment #%d:%s" % (pole1segnum,str(segpole1)))

		#Area2 => Descending phase
		
		idx_pole2   = segpole2.si
		self.log.info("Second poleward boundary at segment #%d:%s" % (pole2segnum,str(segpole2)))
		
		#Start by assuming that the equatorward boundary is
		#simply the end of the segment for which we found the second poleward boundary
		#If we don't find anything better equatorward of that we will stick with it
		idx_equator2= segpole2.ei
		#pick eq2 ignoring very narrow spikes

		for i,seg in enumerate(reversed(self.segments[pole2segnum+1:])):
			if ( seg.twidth >= self.MIN_EQ_SPIKE_DT and seg.area_uncert < np.nanmean(self['total_flux_std'])):
				self.log.info("Found second equatorward boundary at segment #%d:%s" % (len(self.segments)-i,str(seg)))
				idx_equator2 = seg.ei #if meets time requirement, define as boundary
				break

		if idx_equator2 is None:
			self.log.warn('No 2nd equatorward boundary found satisfying delta-t > %s' % (str(self.MIN_EQ_SPIKE_DT)))

		
		return idx_pole1,idx_equator1,idx_pole2,idx_equator2,segpole1,segpole2	

def copy_test_data_and_return_cdffn(destdir):
	#Determine where this module's source file is located
	#to determine where to look for the test data
	src_file_dir = os.path.dirname(os.path.realpath(__file__))
	test_data_dir = os.path.join(src_file_dir,'test_data')
	test_cdffn = 'dmsp-f16_ssj_precipitating-electrons-ions_20100529_v1.1.2.cdf'
	shutil.copy(os.path.join(test_data_dir,test_cdffn),
				destdir)
	return os.path.join(destdir,test_cdffn)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="DMSP SSJ auroral boundary identification")

	parser.add_argument("dmsp_number", 
							type=int,
							help='Process SSJ for F## where this argument is ##',
							default=None)
	parser.add_argument("year", 
							type=int,
							help='Year of date to process',
							default=None)
	parser.add_argument("month",
							type=int,
							help='Month of date to process',
							default=None)
	parser.add_argument("day",
							type=int,
							help='Day of date to process',
							default=None)
	parser.add_argument("--nocsv",
							action='store_true',
							help="Don't write a CSV file of results",
							default=False)
	parser.add_argument("--makeplots",
							action='store_true',
							help="Make plots of successful identifiations",
							default=False)
	parser.add_argument("--plotfailed",
							action='store_true',
							help="Plot unsuccessful identifiations",
							default=False)
	parser.add_argument("--dumpdir",
							help="Root directory for CDF data, plots and CSVs",
							default='/tmp/ssj_auroral_boundary')
	parser.add_argument("--test",
							action='store_true',
							help="Test using the included CDF",
							default=False)

	args = parser.parse_args()
	
	# Create a logging handler if we're running as a script
			
	# create console handler with a lower log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	starttime = datetime.datetime.now().strftime('%H:%M:%S')

	#Formatting for the log (shows colors in the log on the console)
	cformatter = ColoredFormatter(
			"%(white)s%(name)s%(reset)s - %(log_color)s%(levelname)s%(reset)s %(white)s%(message)s",
			datefmt='%m-%d %H:%M:%S',
			reset=True,
			log_colors={
					'DEBUG':    'cyan',
					'INFO':     'green',
					'WARNING':  'yellow',
					'ERROR':    'red',
					'CRITICAL': 'red,bg_white',
			},
			secondary_log_colors={},
			style='%'
	)

	# create formatter and add it to the handlers

	ch.setFormatter(cformatter)

	# add the handlers to the logger
	log.addHandler(ch)

	#Output root directory
	dumpdir_root = args.dumpdir
	
	#Expected CDF file
	cdffn = ('dmsp-f%.2d' % (args.dmsp_number)
	 		+'_ssj_precipitating-electrons-ions_'
			+'%d%.2d%.2d_v1.1.2.cdf' % (args.year,args.month,args.day))

	#Use text of CDF filename as subdirectory
	cdffn = os.path.split(cdffn)[-1]
	cdffn_as_subdir = os.path.splitext(cdffn)[0]
	dump_subdir = cdffn_as_subdir
	
	#All directories
	dumpdirs = {
			'data': os.path.join(dumpdir_root,'data'),
			'image': os.path.join(dumpdir_root,dump_subdir,'img'),
			'csv': os.path.join(dumpdir_root,dump_subdir,'csv')
			}
	try:
		for dirkey,directory in dumpdirs.iteritems():
			if not os.path.exists(directory):
				os.makedirs(directory)
	except:
		print(traceback.format_exc())
		log.fatal('Unable to create output directory %s' % (dumpdir_root))
		print(('Could not create output directory %s.\n' % (dumpdir_root)
				+' You can manually specify where to put data'
				+' using the --dumpdir command line argument'))
		raise
	
	#Use the test CDF file, overriding the other settings if we are running
	#the builtin test run
	if args.test:
		cdffn = copy_test_data_and_return_cdffn(dumpdirs['data'])
		log.warn('--test command line argument overrides '
				 'satellite, year, month, day settings')
	else:
		cdffn = os.path.join(dumpdirs['data'],cdffn)

	log.info('Beginning run of CDF file %s' % (cdffn))
	
	absd = absatday(cdffn,
					imgdir=dumpdirs['image'],
					csvdir=dumpdirs['csv'],
					make_plot=args.makeplots,
					plot_failed=args.plotfailed,
					writecsv=(not args.nocsv))


