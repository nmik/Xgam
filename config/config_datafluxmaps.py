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

import os
from Xgam import X_OUT

OUT_LABEL = '12years_SV_t56'

IN_LABELS_LIST = ['1yr_SV_t56', '2yr_SV_t56',
			      '3yr_SV_t56', '4yr_SV_t56',
			      '5yr_SV_t56', '6yr_SV_t56',
			      '7yr_SV_t56', '8yr_SV_t56',
			      '9yr_SV_t56', '10yr_SV_t56',
			      '11yr_SV_t56', '12yr_SV_t56']

# Activated with --foresub True
FORE_FILES_LIST = os.path.join(X_OUT, 'fits/gll_iem_v06_hp512_list.txt') 
FORE_LABEL = 'glliemv06'

BINNING_LABEL = 'mybins'
MICRO_BINS_FILE = os.path.join(X_OUT, 'ebinning.fits')
MACRO_BINS = [(21,24),(25,32),(33,41),(42,50),(51,66),(67,86)]

IGRB_FILE = os.path.join('config/IGRB_2015.txt') #None if foreground subtraction is not performed
BINCALC = 'CENTER' #CENTER or EDGE, matching dataset selection

MASK_FILE = os.path.join(X_OUT, 'fits/Masksmart_hp512_gp25_src4fgl_testbins_list.txt')  #'None' if no mask is needed
MASK_LABEL = 'SMgp25src4fgl'
