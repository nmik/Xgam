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

"""Produces sky masks (work in progress)
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg

__description__ = 'Produce masks fits files'


"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-ff', '--maskfiles', type=str, required=True, nargs='*',
                    help='Mask files to combine.')
PARSER.add_argument('--outflabel', type=str, required=True,
                    help='Mask files to combine.')                    
PARSER.add_argument('--show', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='if True the mask map is displayed')

def combinemasks(**kwargs):
    """Routine to produce sky maps (limited edition)
    """
    f_ = kwargs['maskfiles']
    logger.info('Reading mask %s...'%f_[0])
    combmask = hp.read_map(f_[0])
    
    for f in f_[1:]:
        m = hp.read_map(f)
        combmask *= m
        
    hp.write_map(os.path.join(X_OUT, 'fits/MaskCombo_%s.fits'%kwargs['outflabel']), combmask)
    logger.info('Created %s'%os.path.join(X_OUT, 'fits/MaskCombo_%s.fits'%kwargs['outflabel']))

    if kwargs['show'] == True:
    
    	hp.mollview(combmask, cmap='bone')
    	plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    combinemasks(**args.__dict__)
