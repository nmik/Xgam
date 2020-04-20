#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#

import os
from Xgam import X_OUT

IN_LABELS_LIST = ['w9w50_SV_t56', 'w51w100_SV_t56']

FORE_FILES_LIST = [os.path.join(X_OUT,'fits/gll_iem_v07_hp512_407.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_528.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_687.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_893.fits'),
			       os.path.join(X_OUT,'fits/gll_iem_v07_hp512_1161.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_1509.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_1961.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_2549.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_3313.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_4305.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_7273.fits'),
                   os.path.join(X_OUT,'fits/gll_iem_v07_hp512_12284.fits')]
FORE_LABEL = 'glliemv07'
OUT_LABEL = '100weeks_SV_t56'

BINNING_LABEL = '2binsTEST'
MICRO_BINS_FILE = os.path.join(X_OUT, 'ebinning.fits')
MACRO_BINS = [(27,32),(33, 38)]

POWER_LAW_INDEX = 2.30
IGRB_FILE = os.path.join(X_OUT, 'data/iso_P8R3_SOURCEVETO_V2_v1.txt') #None if foreground subtraction is not performed
BINCALC = 'CENTER' #CENTER or EDGE, matching dataset selection

MASK_FILE = os.path.join(X_OUT, 'fits/Mask_gp25_4fgl_extsrc.fits')  #'None' if no mask is needed
MASK_LABEL = 'mymask'
