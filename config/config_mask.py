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
    >>> python bin/mkmask.py -c config/config_mask.py --srcmask True --gpmask flat 
    To see the option availabe for bin/mkdataselection.py type:
    >>> python bin/mkmask.py -h
"""

import os
from Xgam import X_CONFIG

OUT_LABEL = 'Mask_weighted1_1000.fits'                                                                                                                     

NSIDE = 128                                                                                                                                                      
SRC_CATALOG = '/Users/mnegro/Documents/_FERMI/CATALOGS/gll_psc_v19_4FGL.fit'
EXTSRC_CATALOG = '/Users/mnegro/Documents/_FERMI/CATALOGS/gll_psc_v19_4FGL.fit'
SRC_MASK_RAD = 1 #[deg]                                                                                                                                               
GP_MASK_LAT = 25. #[deg]

# Only for weighted mask usinf function "mask_src_fluxPSFweighted_1"
PSF_FILE = os.path.join(X_CONFIG, 'fits/psf_SV_t32.fits')
ENERGY = 1000 #[MeV]

# Only for Northen/Southern hemisphere masks
NORTH_LAT = 80 # mask from this value up
SOUTH_LAT = -2 # mask from this value down