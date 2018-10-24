# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import os
import glob

os.environ['DISTUTILS_DEBUG'] = "1"

from setuptools import setup, Extension
from setuptools.command import install as _install

setup(name='ssj_auroral_boundary',
      version = "0.1",
      description = "Figure of Merit Boundary Identification for DMSP SSJ5",
      author = "Liam Kilcommons",
      author_email = 'liam.kilcommons@colorado.edu',
      url = "https://github.com/lkilcommons/ssj_auroral_boundary",
      download_url = "https://github.com/lkilcommons/ssj_auroral_boundary",
      long_description = \
            """
            This package implements the Figure-of-Merit auroral boundary
            identification technique for electron precipitation data from 
            the SSJ Precipitating Electrons and Ions instrument aboard 
            Defense Meteorology Satellite Program (DMSP) spacecraft.
            The algorithm here is tuned for DMSP carrying version 5 of the SSJ
            detector, that is DMSP F16, F17 and F18.
            """,
      install_requires=['numpy','matplotlib','spacepy','colorlog'],
      packages=['ssj_auroral_boundary'],
      package_dir={'ssj_auroral_boundary' : 'ssj_auroral_boundary'},
      package_data={'ssj_auroral_boundary': ['test_data/*']},
      license='LICENSE.txt',
      zip_safe = False,
      classifiers = [
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Programming Language :: Python"
            ],
      )
