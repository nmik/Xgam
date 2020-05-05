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

"""Produces sky masks for the macro-bins, returns the PSF file for the whole
dataset, computes the livetime merging if necessary. Takes a configuration
file similar to the fluxmaps script. More parameters can be defined with
keyword arguments (see --help).
"""

import os
import ast
import sys
import argparse
import numpy as np
import healpy as hp
from astropy.io import fits as pf


from Xgam import X_OUT
from Xgam import X_CONFIG
from Xgam import FT_DATA_FOLDER
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.parsing_ import get_energy_from_fits

__description__ = 'Produce masks fits files.'

"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--configpsf', type=str, required=False,
                    help='PSF configuration file')
PARSER.add_argument('--psffile', type=str, required=False,
                    help='PSF input file', default='default_psf.fits')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='sources mask activated')
#PARSER.add_argument('--gpmask', type=str, choices=['no','shape', 'flat'],
#                    default='no',
#                    help='galactic plain mask (only "flat" available now)')
PARSER.add_argument('--gpcut',type=ast.literal_eval,
                    default=25.0, help='Galactic cut (deg)')
PARSER.add_argument('--nside',type=ast.literal_eval,
                    default=1024, help='NSIDE of output maps')
PARSER.add_argument('--srccat',type=str,
                    default='gll_psc_v20.fit', help='sources catalog in fits folder')
PARSER.add_argument('--srcextcat',type=str,
                    default='gll_psc_v20.fit', help='extented sources catalog in fits folder')
PARSER.add_argument('--outflabel',type=str,
                    default='4FGL_GP25', help='label of output file')
PARSER.add_argument('--typesrcmask', type=ast.literal_eval,
                    choices=[1, 2], default=1,
                    help='type of weighted sources mask')
PARSER.add_argument('--show', type=bool, choices=[True, False],
                    default=False, help='if True the mask map is displayed')
PARSER.add_argument('--overwrite', type=bool, choices=[True, False],
                    default=False, help='if True overwrite existing masks')

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

def mkSmartMask(**kwargs):
    """Routine to produce sky maps, according to defined macro-bins.
    """
    logger.info('SmartMask production...')
    logger.info('Checking PSF file...')
    #get_var_from_file(kwargs['config'])

    PSF_FILE = kwargs['psffile']
    PSF_FILE = os.path.join(FT_DATA_FOLDER, 'fits/' + PSF_FILE)
    if os.path.exists(PSF_FILE):
        logger.info('OK')
    else:
        logger.info('ATT: File not found: %s'%PSF_FILE)
        logger.info('Creating PSF file...')
        assert(kwargs['configpsf'].endswith('.py'))
        get_var_from_file(kwargs['configpsf'])
        try:
            LT_FILE = data.LT_FILE
        except:
            LT_FILE = ''
        if os.path.exists(LT_FILE):
            logger.info('Livetime file found.')
        else:
            try:
                label_lt = data.OUT_LABEL_LT
            except:
                label_lt = 'default'
            lt_dict = data.SUMLT_DICT
            from Xgam.utils.ScienceTools_ import sumLT
            out_gtltsum = sumLT(label_lt,lt_dict)
        from Xgam.utils.ScienceTools_ import gtpsf
        try:
            label_psf = data.OUT_LABEL_PSF
        except:
            label_psf = 'default'
        psf_dict = data.PSF_DICT
        PSF_FILE = gtpsf(label_psf,psf_dict)

    assert(kwargs['config'].endswith('.py'))
    get_var_from_file(kwargs['config'])
    macro_bins = data.MACRO_BINS
    mask_label = data.MASK_LABEL
    in_labels_list = data.IN_LABELS_LIST
    micro_bin_file = data.MICRO_BINS_FILE

    #NOTE: IMPLEMENT TYPE2 ??
    if kwargs['typesrcmask'] == 1:
        from Xgam.utils.mkmask_ import mask_src_fluxPSFweighted_1 as mask_src
    from Xgam.utils.mkmask_ import mask_gp
    from Xgam.utils.parsing_ import get_psf_en_univariatespline

    nside = kwargs['nside']
    out_label = kwargs['outflabel']
    npix = hp.nside2npix(nside)
    src_cat = os.path.join(FT_DATA_FOLDER, 'fits/' + kwargs['srccat'])
    src_ext_cat = os.path.join(FT_DATA_FOLDER, 'fits/' + kwargs['srcextcat'])

    for i, (minb, maxb) in enumerate(macro_bins):
        logger.info('Considering bins from %i to %i...' %(minb, maxb))
        emin, emax, emean = get_energy_from_fits(micro_bin_file, minbinnum=minb, maxbinnum=maxb)
        E_MIN, E_MAX = emin[0], emax[-1]
        E_MEAN = np.sqrt(emax[0]*emin[-1])
        logger.info('Energies %.1f - %.1f [MeV]' %(E_MIN, E_MAX))

        out_name = os.path.join(X_OUT, 'fits/Mask_'+out_label+'_%.1f_MeV.fits'%(E_MIN))
        if os.path.exists(out_name) and not kwargs['overwrite']:
            logger.info('ATT: %s already exists! \n'%out_name)
        else:
            bad_pix = []
            mask = np.ones(npix)
            bad_pix += mask_gp(kwargs['gpcut'], nside)
            psf_spline = get_psf_en_univariatespline(PSF_FILE)
            bad_pix += mask_src(src_cat, src_ext_cat, psf_spline, E_MIN, nside)

            for bpix in np.unique(bad_pix):
                mask[bpix] = 0
            if not os.path.exists(os.path.join(X_OUT, 'fits')):
                os.system('mkdir %s' %os.path.join(X_OUT, 'fits'))
            fsky = 1-(len(np.unique(bad_pix))/float(npix))
            logger.info('fsky = %.3f'%fsky)
            hp.write_map(out_name, mask, coord='G')
            logger.info('Created %s \n' %out_name)

            if kwargs['show'] == True:
                import matplotlib.pyplot as plt
                hp.mollview(mask, cmap='bone')
                plt.show()

        '''
        #OLD IDEA, NOT ELEGANT
        config_temp = os.path.join(FT_DATA_FOLDER, 'conf_mask.py')
        data_mask_dict = {'NSIDE' : kwargs['nside'],
                          'OUT_LABEL' : kwargs['outflabel'],
                          'SRC_CATALOG' : src_cat,
                          'EXTSRC_CATALOG': src_ext_cat,
                          'GP_MASK_LAT' : kwargs['gpcut'],
                          'PSF_FILE' : PSF_FILE,
                          'ENERGY' : np.round(E_MIN,1)}
        mask_dict = {'config': config_temp,
                 'srcmask' : False,
                 'extsrcmask' : False,
                 'gpmask' : 'flat',
                 'srcweightedmask' : True,
                 'srcweightedmask2' : False,
                 'northmask' : False,
                 'southmask' : False
                 }
        print(data_mask_dict)
        print(mask_dict)
        logger.info('Creating temporary config file for mask...')
        createConfigMask(config_temp,data_mask_dict)

        '''
    logger.info('Done!')

#def createConfigMask(config_out,data_mask):
#    f = open(config_out, "w")
#    f.write('from Xgam import X_CONFIG\n\n')
#    for k, v in data_mask.items():
#        if type(v) == str:
#            f.write(str(k) + ' = \'' + str(v) + '\'\n')
#        else:
#            f.write(str(k) + ' = ' + str(v) + '\n')
#    f.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkSmartMask(**args.__dict__)
