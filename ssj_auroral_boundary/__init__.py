# Copyright 2018 SEDA Group at CU Boulder
# Created by:
# Liam Kilcommons
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
"""
ssj_auroral_boundary
--------------------

Figure of Merit boundary identification for DMSP SSJ5

Modules
-------
absatday
abpolarpass
absegment
abscv
files
dmsp_spectrogram

"""
from __future__ import print_function

__version__ = str("0.1.1")

#Prefix for all package loggers
loggername = 'ssj_auroral_boundary'

__all__ = ['absatday', 'abpolarpass', 'absegment', 'abcsv', 'files',
           'dmsp_spectrogram']

# Explicitly import all modules (in addition to defining __all__)
from ssj_auroral_boundary import (absatday,
                                    abpolarpass,
                                    absegment,
                                    abcsv,
                                    files,
                                    dmsp_spectrogram)
