#!/usr/bin/env python
# Copyright 2018 SEDA Group at CU Boulder
# Created by:
# Liam Kilcommons
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # Use the non-GUI backend
from matplotlib import pyplot as pp
import matplotlib.transforms as mtransforms
import logging,sys,datetime,os,argparse,shutil,traceback

from ssj_auroral_boundary import files
from ssj_auroral_boundary.absatday import absatday

from ssj_auroral_boundary import loggername
log = logging.getLogger(loggername)
log.setLevel(logging.DEBUG)

from colorlog import ColoredFormatter #This is from PyPI

def prepare_output_directories(data_rootdir,out_subdir):
    dirs = {
            'data': os.path.join(data_rootdir,'data'),
            'image': os.path.join(data_rootdir,out_subdir,'img'),
            'csv': os.path.join(data_rootdir,out_subdir,'csv')
            }

    for dirkey,directory in dirs.items():
        if not os.path.exists(directory):
            os.makedirs(directory)

    return dirs

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="DMSP SSJ auroral boundary identification")

    parser.add_argument("dmsp_number",
                            type=int,
                            help='Process SSJ for F## where this argument is ##',
                            default=None)
    parser.add_argument("year",
                            type=int,
                            help='Year of date to process',
                            default=None)
    parser.add_argument("month",
                            type=int,
                            help='Month of date to process',
                            default=None)
    parser.add_argument("day",
                            type=int,
                            help='Day of date to process',
                            default=None)
    parser.add_argument("--nocsv",
                            action='store_true',
                            help="Don't write a CSV file of results",
                            default=False)
    parser.add_argument("--makeplots",
                            action='store_true',
                            help="Make plots of successful identifiations",
                            default=False)
    parser.add_argument("--plotfailed",
                            action='store_true',
                            help="Plot unsuccessful identifiations",
                            default=False)
    parser.add_argument("--datarootdir",
                            help="Root directory for CDF data, plots and CSVs",
                            default='/tmp/ssj_auroral_boundary')
    parser.add_argument("--test",
                            action='store_true',
                            help="Test using the included CDF",
                            default=False)
    parser.add_argument("--quiet",
                            action='store_true',
                            help="Suppress most log messages (loglevel=FATAL)")

    args = parser.parse_args()

    # Create a logging handler if we're running as a script

    # create console handler with a lower log level
    ch = logging.StreamHandler()

    loglevel = logging.FATAL if args.quiet else logging.INFO
    ch.setLevel(loglevel)

    starttime = datetime.datetime.now().strftime('%H:%M:%S')

    #Formatting for the log (shows colors in the log on the console)
    cformatter = ColoredFormatter(
            "%(white)s%(name)s%(reset)s - %(log_color)s%(levelname)s%(reset)s %(white)s%(message)s",
            datefmt='%m-%d %H:%M:%S',
            reset=True,
            log_colors={
                    'DEBUG':    'cyan',
                    'INFO':     'green',
                    'WARNING':  'yellow',
                    'ERROR':    'red',
                    'CRITICAL': 'red,bg_white',
            },
            secondary_log_colors={},
            style='%'
    )

    # create formatter and add it to the handlers

    ch.setFormatter(cformatter)

    # add the handlers to the logger
    log.addHandler(ch)

    #Output root directory
    data_rootdir = args.datarootdir

    #Spacecraft and date
    dmsp_number = args.dmsp_number
    year,month,day = args.year,args.month,args.day

    #CSV/Plot subdirectory for this spacecraft
    out_subdir = 'F%.2d_%d%.2d%.2d'%(dmsp_number,year,month,day)

    try:
        dirs = prepare_output_directories(data_rootdir,out_subdir)
    except:
        print(traceback.format_exc())
        log.fatal('Could not create output directories.' % (data_rootdir))
        raise

    if args.test:

        #Use the test CDF file, overriding the satnum and date settings
        test_data_dir,test_cdffn = files.test_cdf_path_and_filename()
        cdffn = os.path.join(test_data_dir,test_cdffn)
        shutil.copy(cdffn,dirs['data'])

        log.warn('--test command line argument overrides '
                 'satellite, year, month, day settings')
    else:

        #Try to download CDF if can't find it in the data directory
        url_and_fn = files.cdf_url_and_filename(dmsp_number,year,month,day)
        cdf_url,cdffn = url_and_fn
        cdffn = os.path.join(dirs['data'],cdffn)
        if not os.path.exists(cdffn):
            log.info('CDF %s does not exist, attempting to download' % (cdffn))
            files.download_cdf_from_noaa(cdf_url,cdffn)

    log.info('Beginning run of CDF file %s' % (cdffn))

    absd = absatday(cdffn,
                    imgdir=dirs['image'],
                    csvdir=dirs['csv'],
                    make_plot=args.makeplots,
                    plot_failed=args.plotfailed,
                    writecsv=(not args.nocsv))


