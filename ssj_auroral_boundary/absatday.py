# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import numpy as np
import logging,sys,datetime,os,traceback
from geospacepy import special_datetime
from spacepy import pycdf

from ssj_auroral_boundary import loggername
from ssj_auroral_boundary.abpolarpass import abpolarpass
from ssj_auroral_boundary.abcsv import abcsv

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
		self.csv = abcsv(csvdir,csvfile,cdffile,writecsv=self.writecsv)

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
		   