.. ssj_auroral_boundary documentation master file, created by
   sphinx-quickstart on Wed Oct 24 15:39:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ssj_auroral_boundary's documentation!
================================================

.. toctree::
   :maxdepth: 2

   quickstart
   api

This package implements an algorithm for identifying when and where
the Defense Meteorology Satellite Program (DMSP) spacecraft encounter the
auroral oval.

The algorithm produces time and location (in magnetic coordinates)
of the equatorward and polarward edges of the aurora based on
changes in the energy and intensity of precipitating electrons 
measured by the Special Sensor for Electrons and 
Ions Version 5 (SSJ5) instrument.

These auroral boundaries are in-situ, not remotely sensed, meaning that
they only provide information about the auroral oval at the spacecraft's
location.

The algorithm is further described in:

	Kilcommons, L. M., Redmon, R. J., & Knipp, D. J. (2017). 
	A New DMSP Magnetometer & Auroral Boundary Dataset and Estimates of Field 
	Aligned Currents in Dynamic Auroral Boundary Coordinates. 
	Journal of Geophysical Research: Space Physics, 2016JA023342. 
	<https://doi.org/10.1002/2016JA023342>

And is based on an earlier boundary identification system described in:

	Redmon, R. J., Peterson, W. K., Andersson, L., Kihn, E. A., 
	Denig, W. F., Hairston, M., & Coley, R. (2010). 
	Vertical thermal O+ flows at 850 km in dynamic auroral boundary coordinates.
	Journal of Geophysical Research: Space Physics, 115(A11), A00J08. 
	<https://doi.org/10.1029/2010JA015589>


