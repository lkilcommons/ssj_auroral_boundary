.. ssj_auroral_boundary documentation master file, created by
   sphinx-quickstart on Wed Oct 24 15:39:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ssj_auroral_boundary's documentation!
================================================

This package implements an algorithm for identifying when and where
the Defense Meteorology Satellite Program (DMSP) spacecraft encounter the
auroral oval.

.. toctree::
   :maxdepth: 2

   quickstart
   api

Data Source
-----------

DMSP Special Sensor 'J' (SSJ) Precipitating Electrons and Ions Instrument


Data Format
-----------

NASA Common Data Format (CDF)


Data Used
---------

* Precipitating electron energy flux (integrated across 9 SSJ channels ~1-30 keV)
* Root sum of squared channel uncertainty (poisson / digitization, same channels)

.. rubric:: Table of SSJ Channel Center Energies [keV]

===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
  0     1     2     3     4     5     6     7     8     9 
===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
30.18 20.62 14.04 9.58  6.50  4.42  3.05  2.06  1.41  .992 	
===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

Default Output
--------------

* Comma Seperated Value (CSV) text file
* One file per spacecraft, per day
* Universal Time, Magnetic Latitude, Magnetic Local Time for each boundary
* :download:`Sample CSV <_static/sample_output.csv>`

Optional Graphics
-----------------

* One graphic per polar pass (2 auroral crossings)
* :download:`Sample PNG <_static/sample_plot.png>`


Further Reading
---------------

* Current algorithm described in:

	Kilcommons, L. M., Redmon, R. J., & Knipp, D. J. (2017). 
	A New DMSP Magnetometer & Auroral Boundary Dataset and Estimates of Field 
	Aligned Currents in Dynamic Auroral Boundary Coordinates. 
	Journal of Geophysical Research: Space Physics, 2016JA023342. 
	https://doi.org/10.1002/2016JA023342

* Figure-of-Merit boundaries first described in:

	Redmon, R. J., Peterson, W. K., Andersson, L., Kihn, E. A., 
	Denig, W. F., Hairston, M., & Coley, R. (2010). 
	Vertical thermal O+ flows at 850 km in dynamic auroral boundary coordinates.
	Journal of Geophysical Research: Space Physics, 115(A11), A00J08. 
	https://doi.org/10.1029/2010JA015589


* DMSP SSJ data version used described in:

	Redmon, R. J., Denig, W. F., Kilcommons, L. M., & Knipp, D. J. (2017). 
	New DMSP Database of Precipitating Auroral Electrons and Ions. 
	Journal of Geophysical Research: Space Physics, 2016JA023339. 
	https://doi.org/10.1002/2016JA023339

