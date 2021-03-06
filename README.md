# ssj_auroral_boundary
Identify boundaries of the aurora with Defense Meteorology Satellite Program (DMSP) electron precipitation.

[![Documentation Status](https://readthedocs.org/projects/ssj-auroral-boundary/badge/?version=latest)](https://ssj-auroral-boundary.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/154356043.svg)](https://zenodo.org/badge/latestdoi/154356043)

## How to cite

* Cite the Zenodo DOI (above) for the code
* Cite the paper (below) for the algorithm

## Algorithm

This code implements the auroral boundary identification technique described in:

Kilcommons, L. M., Redmon, R. J., & Knipp, D. J. (2017). A New DMSP Magnetometer & Auroral Boundary Dataset and Estimates of Field Aligned Currents in Dynamic Auroral Boundary Coordinates. Journal of Geophysical Research: Space Physics, 2016JA023342. https://doi.org/10.1002/2016JA023342

## I really just want some boundary data

The auroral boundaries (as daily CSV files) for F16,F17 and F18 2010 - 2014 are on Zenodo. Please cite the following DOI if you use the dataset in a publication, as well as citing the paper above.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3373812.svg)](https://doi.org/10.5281/zenodo.3373812)

## Installation

## Tested Operating Systems
* Ubuntu Linux (14.04 and 18.04)
* MacOS (Sierra, High Sierra)

## Tested Python Versions
* 2.7 (Anaconda2)
* 3.6 (Anaconda)
* 3.8 (Anaconda)

### Other Software Required
There are a few tools used by this library that can't be installed as part of
the main install script (setup.py).

1. [NASA CDF library](https://cdf.gsfc.nasa.gov/html/sw_and_docs.html)
This software expects the DMSP data it uses to be packaged as 
NASA Common Data Format (CDF) files. This library is required to read those
files. An additional Python package (Spacepy) provides the python interface.

2. [Geospacepy](https://github.com/lkilcommons/geospacepy-lite)
This is a small python package which provides plotting routines and datetime
handling.

## Clone the github repo
```
git clone https://github.com/lkilcommons/ssj_auroral_boundary.git
```

## Run the installer
```
cd ssj_auroral_boundary
python setup.py install
```

## HELP!
If something went wrong or you've found a bug or incompatability please file
an issue (above).
