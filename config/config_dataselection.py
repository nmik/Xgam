#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#

""" Example of configuration file for bin/mkmask.py.

    To run the analysis just do:
    >>> python bin/mkmask.py -c config/config_mask.py
    To see the option availabe for bin/mkdataselection.py type:
    >>> python bin/mkdataselection.py -h
"""


import os
import numpy as np
from Xgam import FT_DATA_FOLDER
from Xgam.utils.ScienceTools_ import gtEbindef


OUT_LABEL = 'w9w50_SV_t32'

#start and stop week number to be selected in FT1 and FT2 files
START_WEEK = 9
STOP_WEEK = 50

# minumum and maximum energy and energy bin number to generate fits file
E_MIN = 1000
E_MAX = 1000000
E_NBINS = 2
EBINNING_ARRAY = np.logspace(np.log10(E_MIN), np.log10(E_MAX), E_NBINS+1)
EBINNING_FILE = gtEbindef(EBINNING_ARRAY, file_name='ebinning.txt')
SC_FILE = os.path.join(FT_DATA_FOLDER, 'spacecraft', 'lat_spacecraft_merged.fits')

# other settings
ZMAX = 90
EVCLASS = 2048
EVTYPE = 32
IRFS = 'P8R3_SOURCEVETO_V2'
HPX_MAP_ORDER = 5
FILTER_CUT='DATA_QUAL==1&&LAT_CONFIG==1&&LAT_MODE==5&&IN_SAA!=T'+\
               '&&((ABS(ROCK_ANGLE)<52))'

GTSELECT_DICT = {'infile': 'DEFAULT',
                 'emin': E_MIN,
                 'emax': E_MAX,
                 'zmax': ZMAX,
                 'evclass': EVCLASS,
                 'evtype': EVTYPE,
                 'outfile': 'DEFAULT',
                 'chatter': 4,
                 'clobber': 'no',
                 'tmin' : 'INDEF',
                 'tmax': 'INDEF'}

GTMKTIME_DICT = {'evfile': 'DEFAULT',
                 'scfile': SC_FILE,
                 'filter': FILTER_CUT,
                 'roicut': 'no',
                 'outfile': 'DEFAULT',
                 'clobber': 'no'}

GTBIN_DICT = {'evfile': 'DEFAULT',
              'algorithm': 'HEALPIX',
              'scfile': SC_FILE,
              'hpx_ordering_scheme': 'RING',
              'hpx_order': HPX_MAP_ORDER,
              'coordsys': 'GAL',
              'hpx_ebin': 'yes',
              'ebinalg': 'FILE',
              'ebinfile': EBINNING_FILE,
              'outfile': 'DEFAULT',
              'clobber': 'no'}

GTLTCUBE_DICT = {'evfile': 'DEFAULT',
                 'scfile':  SC_FILE,
                 'zmax': ZMAX,
                 'dcostheta': 0.025,
                 'binsz': 1,
                 'outfile': 'DEFAULT',
                 'chatter': 4,
                 'clobber': 'no'}

GTEXPCUBE2_DICT = {'infile': 'DEFAULT',
                   'cmap': 'DEFAULT',
                   'hpx_ordering_scheme': 'RING',
                   'hpx_order': HPX_MAP_ORDER,
                   'irfs': IRFS,
                   'evtype': EVTYPE,
                   'outfile': 'DEFAULT',
                   'ebinalg': 'FILE',
              	   'ebinfile': EBINNING_FILE,
                   'bincalc' : 'CENTER',
                   'clobber': 'no'}
