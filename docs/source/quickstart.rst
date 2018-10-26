Quickstart Guide
================

Checking your installation
--------------------------

The ssj_auroral_boundary code uses DMSP SSJ data packaged in NASA CDF files.
Before beginning to generate boundaries, ensure that the NASA CDF library is
installed on your computer.

The CDF package can be downloaded from the
`NASA CDF Website <https://cdf.gsfc.nasa.gov/html/sw_and_docs.html>`_

Unfortunately NASA does not provide Python CDF support,
so this code relies on CDF-reading functionality provided by Spacepy.

A sample CDF file is included with the package. Running this code,

.. testcode::

	import os
	from ssj_auroral_boundary import files
	from spacepy import pycdf
	path,filename = files.test_cdf_path_and_filename()
	with pycdf.CDF(os.path.join(path,filename)) as cdf:
		print(cdf.attrs['Generated_by'])

Should generate,

.. testoutput::

	Rob Redmon

Ways to generate boundaries
---------------------------


Command Line Interface
++++++++++++++++++++++

Boundaries are most easily generated using the command line interface.
After installing the package, a script called find_ssj_boundaries.py is added
to your python path.

You can run the boundary identification algorithm using the test CDF file using
the terminal:

.. code:: bash
	
	find_ssj_boundaries.py 16 2012 1 1 --test

.. note::

	Notice that you don't need to put a :code:`python` before
	:code:`find_ssj_boundaries.py`. While it does have the .py extension
	this script is just a command as far as your terminal is concerned.
	It executes using whatever python interpreter would be the default 
	if you typed :code:`python` into your terminal.

This should produce a large number of logging messages on your terminal:

.. code:: bash

	ssj_auroral_boundary - WARNING --test command line argument overrides satellite, year, month, day settings
	ssj_auroral_boundary - INFO Beginning run of CDF file dmsp-f16_ssj_precipitating-electrons-ions_20100529_v1.1.2.cdf
	ssj_auroral_boundary.absatday - INFO Satellite number determined to be 16
	ssj_auroral_boundary.absatday - WARNING Unable to find APEX latitude or local time variables in CDF file. Falling back to AACGM magnetic coordinates
	ssj_auroral_boundary.absatday - INFO Entered Southern Hemisphere: ind:316,lat:-20.126
	ssj_auroral_boundary.absatday - INFO Entered Northern Hemisphere: ind:3309,lat:20.176
	ssj_auroral_boundary.absatday - INFO Entered Southern Hemisphere: ind:6435,lat:-20.126
	ssj_auroral_boundary.absatday - INFO Entered Northern Hemisphere: ind:9514,lat:20.174
	ssj_auroral_boundary.absatday - INFO Entered Southern Hemisphere: ind:12564,lat:-20.126
	ssj_auroral_boundary.absatday - INFO Entered Northern Hemisphere: ind:15706,lat:20.172
	ssj_auroral_boundary.absatday - INFO Entered Southern Hemisphere: ind:18640,lat:-20.125
	ssj_auroral_boundary.absatday - INFO Entered Northern Hemisphere: ind:21885,lat:20.170
	ssj_auroral_boundary.absatday - INFO Entered Southern Hemisphere: ind:24699,lat:-20.122

To suppress these messages you can add the :code:`--quiet` flag to the above command

By default the script creates CSV files in the /tmp directory. 
(the :code:`--datarootdir /path/to/some/place` option changes this)

.. code:: bash

	cd /tmp/ssj_auroral_boundary/F16_20100529/csv
	head dmsp-f16_ssj_precipitating-electrons-ions_20100529_v1.1.2_boundaries.csv

.. code:: bash
	
	# DMSP SSJ Auroral Boundary Identification (dmsp-f16_ssj_precipitating-electrons-ions_20100529_v1.1.2.cdf)
	# Generated on Fri Oct 26 13:52:58 2018
	# Glossary:
	# EQ1: First equator-side auroral boundary 
	# PO1: First pole-side auroral boundary 
	# PO2: Second pole-side auroral boundary 
	# EQ2: Second equator-side auroral boundary 
	# FOM: Figure of merit / Quality ( < 1.8 usually questionable)
	# ind: row of DMSP SSJ CDF file corresponding to boundary
	# hemisphere : -1 for southern hemisphere +1 for northern

You can also produce plots (saved to same directory) by adding the 
:code:`--makeplots` option.


Directly with Python
++++++++++++++++++++

The API guide explains more about how the library works internally,
but this is the minimum code you need to get one set of boundaries' 
time, magnetic latitude and magnetic local time.

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

.. testoutput::
	
	('Sat May 29 02:58:07 2010', 65.541316005828321, 16.88996581629161)
	('Sat May 29 03:00:16 2010', 70.908033561060606, 16.147218793047333)
	('Sat May 29 03:08:38 2010', 72.781246453140852, 9.9421193209708267)
	('Sat May 29 03:11:36 2010', 65.473878712612631, 8.7777971568744455)


