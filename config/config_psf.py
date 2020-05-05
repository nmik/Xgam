#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


""" Example of configuration file for PSF building.

    This should be coupled with mksmartmask.py

    MAYBE MERGE WITH DATA SELECION??
"""

import os
from Xgam import X_CONFIG

OUT_LABEL_LT = 'SV_evt56'
OUT_LABEL_PSF = 'SV_evt56'

LT_LIST = '/archive/home/sammazza/fermi_data/lt_list.txt'

E_MIN = 1000
E_MAX = 10000
EVTYPE = 56
IRFS = 'P8R3_SOURCEVETO_V2'


SUMLT_DICT = {'infile1' : LT_LIST,
              'outfile' : 'DEFAULT',
              'chatter': 4,
              'clobber': 'no'}

PSF_DICT = {'expcube' : 'DEFAULT',
            'outfile' : 'DEFAULT',
            'irfs' : IRFS,
            'evtype' : EVTYPE,
            'ra' : 45.0,
            'dec' : 45.0,
            'emin' : E_MIN,
            'emax' : E_MAX,
            'nenergies' : 500,
            'thetamax' : 30.0,
            'ntheta' : 100,
            'chatter': 4,
            'clobber': 'no'}
