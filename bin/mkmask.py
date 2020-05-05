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
import sys
import argparse
import numpy as np
import healpy as hp
from astropy.io import fits as pf



from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg

__description__ = 'Produce masks fits files'


"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='sources mask activated')
PARSER.add_argument('--extsrcmask', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='extended sources mask activated')
PARSER.add_argument('--gpmask', type=str, choices=['no','shape', 'flat'],
                    default='no',
                    help='galactic plain mask (only "flat" available now)')
PARSER.add_argument('--srcweightedmask', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='weighted sources mask activated')
PARSER.add_argument('--srcweightedmask2', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='weighted sources mask activated (second version)')
PARSER.add_argument('--northmask', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='Northen hemisphere mask activated')
PARSER.add_argument('--southmask', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='Southern hemisphere mask activated')
PARSER.add_argument('--show', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='if True the mask map is displayed')

if (sys.version_info > (3, 0)):
	from importlib.machinery import SourceFileLoader
	def get_var_from_file(filename):
		f = open(filename)
		global data
		data = SourceFileLoader('data', filename).load_module()
		f.close()
else:
	import imp
	def get_var_from_file(filename):
		f = open(filename)
		global data
		data = imp.load_source('data', '', f)
		f.close()

def mkMask(**kwargs):
    """Routine to produce sky maps (limited edition)
    """
    logger.info('Starting mask production...')
    get_var_from_file(kwargs['config'])
    bad_pix = []
    nside = data.NSIDE
    out_label = data.OUT_LABEL
    npix = hp.nside2npix(nside)
    mask = np.ones(npix)
    if kwargs['srcmask'] == True:
        from Xgam.utils.mkmask_ import mask_src
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        bad_pix += mask_src(cat_file, src_mask_rad, nside)
    if kwargs['extsrcmask'] == True:
        from Xgam.utils.mkmask_ import mask_extsrc
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.EXTSRC_CATALOG
        bad_pix += mask_extsrc(cat_file, src_mask_rad, nside)
    if kwargs['gpmask'] == 'flat':
        from Xgam.utils.mkmask_ import mask_gp
        gp_mask_lat = data.GP_MASK_LAT
        bad_pix += mask_gp(gp_mask_lat, nside)
    if kwargs['srcweightedmask'] == True:
    	from Xgam.utils.mkmask_ import mask_src_fluxPSFweighted_1
    	from Xgam.utils.parsing_ import get_psf_en_univariatespline
    	cat_file = data.SRC_CATALOG
    	ext_cat_file = data.EXTSRC_CATALOG
    	psf_f = data.PSF_FILE
    	energy = data.ENERGY
    	psf_spline = get_psf_en_univariatespline(psf_f)
    	bad_pix += mask_src_fluxPSFweighted_1(cat_file, ext_cat_file, psf_spline, energy, nside)
    if kwargs['srcweightedmask2'] == True:
    	from Xgam.utils.mkmask_ import mask_src_fluxPSFweighted_2
    	from Xgam.utils.parsing_ import get_psf_en_univariatespline
    	cat_file = data.SRC_CATALOG
    	ext_cat_file = data.EXTSRC_CATALOG
    	psf_f = data.PSF_FILE
    	e_min = data.E_MIN
    	e_max = data.E_MAX
    	psf_spline = get_psf_en_univariatespline(psf_f)
    	bad_pix += mask_src_fluxPSFweighted_2(cat_file, ext_cat_file, psf_spline, e_min, e_max, nside)
    if kwargs['northmask'] == True:
        from Xgam.utils.mkmask_ import mask_north
        north_lat = data.NORTH_LAT
        bad_pix += mask_north(north_lat, nside)
    if kwargs['southmask'] == True:
        from Xgam.utils.mkmask_ import mask_south
        south_lat = data.SOUTH_LAT
        bad_pix += mask_south(south_lat, nside)
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0
    if not os.path.exists(os.path.join(X_OUT, 'fits')):
    	os.system('mkdir %s' %os.path.join(X_OUT, 'fits'))
    out_name = os.path.join(X_OUT, 'fits/'+out_label+'.fits')
    fsky = 1-(len(np.unique(bad_pix))/float(npix))
    logger.info('fsky = %.3f'%fsky)
    hp.write_map(out_name, mask, coord='G')
    logger.info('Created %s' %out_name)

    if kwargs['show'] == True:
    	import matplotlib.pyplot as plt
    	hp.mollview(mask, cmap='bone')
    	plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkMask(**args.__dict__)
