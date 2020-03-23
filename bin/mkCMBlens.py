
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

"""Produces CMB lensing maps (work in progress)
"""

import os
import argparse
import numpy as np
import healpy as hp
from astropy.io import fits as pf
import cmath

from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg

__description__ = 'Creator of CMB lensing maps from alm measured by Planck Collaboration.'

"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infile', type=str, required=True,
                    help='input fits file of CMB lensing alm')
PARSER.add_argument('-o', '--outlabel', type=str, required=False,
                    help='label of CMB output map', default='CMB_lens')
PARSER.add_argument('--nsideout', type=int, default=512,
                    help='NSIDE of the output HEALPix maps')
PARSER.add_argument('-s', '--show', type=bool, required=False, default=False,
                    help='show produced map')
PARSER.add_argument('--overwrite', type=bool, required=False, default=False,
                    help='overwrite existing map')

def CMB_lens_map_creator(**kwargs):
    """Viewer interface for healpix maps
    """
    input_file = kwargs['infile']
    nside_out = kwargs['nsideout']
    out_label = kwargs['outlabel']
    overwrite = kwargs['overwrite']
    if not os.path.exists(input_file):
        abort("File %s not found!"%input_file)

    logger.info('Reading CMB data...')
    CMB_file = pf.open(input_file) #mean field data from Minimum Variance
    CMB_data = CMB_file[1].data
    logger.info('Filling alm array...')
    real = mf['real']
    imag = mf['imag']
    alm = real + 1j *imag
    logger.info('Generating map with NSIDE %i'%nside_out)
    conv_map = hp.alm2map(alm,nside_out)

    if not os.path.exists(os.path.join(X_OUT, 'fits')):
    	os.system('mkdir %s' %os.path.join(X_OUT, 'fits'))
    out_name = os.path.join(X_OUT, 'fits/'+out_label+'_hp%i.fits'%nside_out)
    hp.write_map(out_name, conv_map, coord='G', overwrite=overwrite)
    logger.info('Created %s'%out_name)

    if kwargs['show'] ==  True:
        import matplotlib as mpl
        mpl.use('Agg')
    	import matplotlib.pyplot as plt
    	mm = hp.ud_grade(conv_map, nside_out=64)
        hp.mollview(mm)
    	plt.show()
    frmaps.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    CMB_lens_map_creator(**args.__dict__)
