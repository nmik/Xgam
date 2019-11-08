#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#

"""Xgam: Framework for Gamma-ray X-correlation Analysis
"""

import os

PACKAGE_NAME = 'Xgam'

"""Basic folder structure of the package.
"""
X_ROOT = os.path.abspath(os.path.dirname(__file__))
X_BIN = os.path.join(X_ROOT, 'bin')
X_CONFIG = os.path.join(X_ROOT, 'config')
X_UTILS = os.path.join(X_ROOT, 'utils')


""" This is where we put the actual (FT1 and FT2) data sets.  
"""

from Xgam.utils.logging_ import logger
try:
    FT_DATA_FOLDER = os.environ['P8_DATA']
    logger.info('Base data folder set to $P8_DATA = %s...' % FT_DATA_FOLDER)
except KeyError:
    FT_DATA_FOLDER = '/Users/mnegro/data/Fermi-LAT'
    logger.info('$P8_DATA not set, base data folder set to %s...' %\
                FT_DATA_FOLDER)

""" This is the output directory.
"""
try:
    X_OUT = os.environ['X_OUT']
    X_OUT_FIG = os.environ['X_OUT_FIG']
except:
    X_OUT = os.path.join(X_ROOT, 'output')
    X_OUT_FIG = os.path.join(X_ROOT, 'output/figures')

if __name__ == '__main__':
	from Xgam.utils.logging_ import startmsg
	startmsg()
	print('X_ROOT: %s' % X_ROOT)
