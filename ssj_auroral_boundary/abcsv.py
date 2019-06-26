# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import os
import datetime as dt
import logging

class abcsv(object):
    """Class for writing satellite day of boundary identifications to CSV file
    """
    def __init__(self, csvdir, csvfn, ssjcdffn, csvvars=['mlat', 'mlt'],
                 writecsv=True):
        """Write output CSV file

        Parameters
        ----------
        csvdir : str
            Output directory (must exist)
        csvfn : str
            Output CSV filename (must exist)
        ssjcdffn : str
            DMSP SSJ CDF filename we are create boundaries for (used in header)
        csvvars : list
            List of variables to write on each line, for each boundary
            (i.e. if mlat, 4 columns will be added per row:
            (mlat of equatorward crossing 1, mlat of polar crossing 1, 
            mlat of polar crossing 2, mlat of equatorward crossing 2)
            Options include mlat (magnetic apex latitude, unless log warning
            notes that AACGM V1 was used), mlt (magnetic apex local time, unless
            log warning notes that AACGM V1 was used), diff_flux, and any valid
            CDF variable names (e.g., SC_GEOCENTRIC_LAT, SC_GEOCENTRIC_LON).
        writecsv : bool
            Calls to add_auroral_boundary_to_csv actually modify a file or not
            Default (True)

        """
        self.csvfn = os.path.join(csvdir,csvfn)
        self.csvvars = csvvars
        self.ssjcdffn = ssjcdffn
        self.writecsv = writecsv
        if self.writecsv:
            self.write_header()

    def write_header(self):
        #Write CSV file header
        
        with open(self.csvfn, 'w') as f:
            #Description
            cdffn = os.path.split(self.ssjcdffn)[-1]
            header_lines = [
                'DMSP SSJ Auroral Boundary Identification (%s)'%(cdffn),
                'Generated on %s' % (dt.datetime.now().strftime('%c')),
                'Glossary:',
                'EQ1: First equator-side auroral boundary ',
                'PO1: First pole-side auroral boundary ',
                'PO2: Second pole-side auroral boundary ',
                'EQ2: Second equator-side auroral boundary ',
                'FOM: Figure of merit / Quality ( < 1.8 usually questionable)',
                'ind: row of DMSP SSJ CDF file corresponding to boundary',
                'hemisphere : -1 for southern hemisphere +1 for northern']
            for line in header_lines:
                f.write('# '+line+'\n')

            #Always added columns
            colnames = ('UTSecond Pass Start,'
                        +'UTSecond Pass End,'
                        +'hemisphere,'
                        +'UTSec EQ1,'
                        +'UTSec PO1,'
                        +'UTSec PO2,'
                        +'UTSec EQ2,'
                        +'FOM,')
            #Optional columns defined in self.csvvars
            for var in self.csvvars:
                for bnd in ['EQ1','PO1','PO2','EQ2']:
                    colnames+="%s %s," % (var,bnd) 
            colnames = colnames[:-1] # remove leftover comma
            colnames += '\n' # add newline
            f.write(colnames)

    def add_auroral_boundary_to_csv(self,abpp):
        """Add a boundary idenfication to CSV

        Parameters
        ----------
        abpp : abpolarpass
            Polar pass object (successful boundary identification)
        
        """
        if self.writecsv:
            hemicode = 1 if abpp.hemi=='N' else -1
            with open(self.csvfn,'a') as f:
                line = ("%d,%d," % (abpp['uts'][0],abpp['uts'][-1])
                        +"%d," % (hemicode)
                        +"%d," % (int(abpp['uts'][abpp.idx_equator1]))
                        +"%d," % (int(abpp['uts'][abpp.idx_pole1]))
                        +"%d," % (int(abpp['uts'][abpp.idx_pole2])) 
                        +"%d," % (int(abpp['uts'][abpp.idx_equator2]))
                        +"%.3f" % (abpp.max_fom))
                for var in self.csvvars:
                    line += ",%.3f,%.3f,%.3f,%.3f" % ( \
                                                abpp[var][abpp.idx_equator1], \
                                                abpp[var][abpp.idx_pole1], \
                                                abpp[var][abpp.idx_pole2], \
                                                abpp[var][abpp.idx_equator2])
                line += "\n"
                f.write(line)
