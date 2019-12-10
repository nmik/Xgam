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


""" Example of configuration file for bin/mkmask.py.

    To run the analysis just do: 
    >>> python bin/mkmask.py -c config/config_mask.py --srcmask True --gpmask flat 
    To see the option availabe for bin/mkdataselection.py type:
    >>> python bin/mkmask.py -h
"""

import os
from GRATools import GRATOOLS_CONFIG

OUT_LABEL = 'Mask_gp25_src4FGL1'#bat'                                                                                                                                       

NSIDE = 128#512                                                                                                                                                       
SRC_CATALOG =  '/Users/mnegro/Documents/_FERMI/CATALOGS/gll_psc_v19_4FGL.fit'
EXTSRC_CATALOG = '/Users/mnegro/Documents/_FERMI/CATALOGS/gll_psc_v19_4FGL.fit'
SRC_MASK_RAD = 1 #[deg]                                                                                                                                               
GP_MASK_LAT = 25.