API Reference
=============

The class structure of this module is based on a parent-child approach. Each
class in this list contains multiple children of the class below

* absatday.absatday - top level class, one spacecraft-day of SSJ data
* abpolarpass.abpolarpass - one half of a DMSP orbit, a crossing of N or S pole
* absegment.absegent - a possible auroral zone, a single span of strong flux

Each successive class is a window to a subset of the parent class's 
DMSP SSJ data. To keep accessing data simple and consistant across all of the
classes, only the top-level class (absatday) stores and names data variables
(time, spacecraft location, and fluxes) The child classes all slice their 
parent's data as needed instead of storing it themselves.

This functionality is implemented with a custom __getitem__ method 
which can be used to fetch the data associated with a specific instance for
any of the 3 classes.

Here is an example:

.. testcode::

	import os
	from ssj_auroral_boundary import files
	path,filename = files.test_cdf_path_and_filename()
	cdffn = os.path.join(path,filename)

	outdir = '/tmp'

	from ssj_auroral_boundary.absatday import absatday

	#Find boundaries
	absd = absatday(cdffn,
	                csvdir=outdir,
	                imgdir=outdir,make_plot=False)

	#Find a successful boundary identification
	good_abpp = None
	for abpp in absd.polarpasses:
	    if abpp.failure_reason is None:
	        good_abpp = abpp
	        break

	#Show how same interface accesses Magnetic Latitude
	day_mlat = absd['mlat']
	pass_mlat = good_abpp['mlat']
	print(day_mlat.shape,pass_mlat.shape)

.. testoutput::
	
	((86400,), (3050,))

Boundaries are represented as indicies into abpolarpass (see abpolarpass below).
You can use this __getitem__ syntax combined with these index variables to
get whatever data value (time,location,flux) you want from among the class
properties in ssj_auroral_boundary.absatday.absatday.

For example following on the previous code snippet:

.. code::

	#Show how to get the boundaries from a successful abpolarpass
	boundary_indicies = [good_abpp.idx_equator1,
	                    good_abpp.idx_pole1,
	                    good_abpp.idx_pole2,
	                    good_abpp.idx_equator2] 

	for idx in boundary_indicies:
	    btime = good_abpp['time'][idx] #datetime
	    bmlat = good_abpp['mlat'][idx] #AACGM magnetic latitude
	    bmlt = good_abpp['mlt'][idx] #AACGM magnetic local time
	    print(btime.strftime('%c'),bmlat,bmlt)

.. code::

	('Sat May 29 02:58:07 2010', 65.541316005828321, 16.88996581629161)
	('Sat May 29 03:00:16 2010', 70.908033561060606, 16.147218793047333)
	('Sat May 29 03:08:38 2010', 72.781246453140852, 9.9421193209708267)
	('Sat May 29 03:11:36 2010', 65.473878712612631, 8.7777971568744455)

	
.. module:: ssj_auroral_boundary
	
.. autoclass:: ssj_auroral_boundary.absatday.absatday

.. autoclass:: ssj_auroral_boundary.abpolarpass.abpolarpass
