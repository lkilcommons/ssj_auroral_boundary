About this Code
===============

This software implements the algorithm described in:

Kilcommons, L. M., Redmon, R. J., & Knipp, D. J. (2017). 
A New DMSP Magnetometer & Auroral Boundary Dataset and Estimates of Field 
Aligned Currents in Dynamic Auroral Boundary Coordinates. 
Journal of Geophysical Research: Space Physics, 2016JA023342. 
https://doi.org/10.1002/2016JA023342

Output
------

* Comma Seperated Value (CSV) text file
* One file per spacecraft, per day
* Universal Time, Magnetic Latitude, Magnetic Local Time for each boundary
* CSV output variables are customizable (see quickstart guide)
* :download:`Sample CSV <_static/sample_output.csv>`

Graphics
--------

* One graphic per polar pass (2 auroral crossings)
* :download:`Sample PNG <_static/sample_plot.png>`

Data Details
------------

Source
++++++

DMSP Special Sensor 'J' (SSJ) Precipitating Electrons and Ions Instrument
Daily NASA CDF file used by this code are automatically downloaded from an 
archive at NOAA (if possible). Data can also be accesed manually at 
NASA CDAWeb.

Format
++++++

NASA Common Data Format (CDF)


Processing
++++++++++

* Precipitating electron energy flux (integrated across 9 SSJ channels ~1-30 keV)
* Root sum of squared channel uncertainty (poisson / digitization, same channels)

.. rubric:: Table of SSJ Channel Center Energies [keV]

===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
  0     1     2     3     4     5     6     7     8     9 
===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
30.18 20.62 14.04 9.58  6.50  4.42  3.05  2.06  1.41  .992 	
===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 


Further Reading
---------------

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
