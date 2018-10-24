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
	from ssj_auroral_boundary import ssj_auroral_boundary
	from spacepy import pycdf
	path,filename = ssj_auroral_boundary.test_cdf_path_and_filename()
	with pycdf.CDF(os.path.join(path,filename)) as cdf:
		print(cdf.attrs['Generated_by'])

Should generate,

.. testoutput::

	Rob Redmon

Ways to generate boundaries
---------------------------

Boundaries are most easily generated using the command line interface.

To run the boundary identification algorithm using the test CDF file.

