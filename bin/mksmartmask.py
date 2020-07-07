#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Simone Ammazzalorso, University of Torino.                            #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
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


from Xgam import X_OUT, FT_DATA_FOLDER 
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
PARSER.add_argument('--irfs', type=str, required=False,
                    help='The Fermi-LAT IRFs to consider')
PARSER.add_argument('--evtype', type=int, required=False,
                    help='The event type to consider')
PARSER.add_argument('--psffile', type=str, required=False,
                    help='PSF input file (if given --irfs and --evtype are discarded)', 
                    default=None)
PARSER.add_argument('--ltfile', type=str, required=False,
                    help='livetimes .fits or .txt with list of .fits', default='')
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
                    help='1=... ; 2=... ;')
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
    assert(kwargs['config'].endswith('.py'))
    get_var_from_file(kwargs['config'])
    macro_bins = data.MACRO_BINS
    mask_label = data.MASK_LABEL
    in_labels_list = data.IN_LABELS_LIST
    micro_bin_file = data.MICRO_BINS_FILE
    emin, emax, emean = get_energy_from_fits(micro_bin_file, \
                        minbinnum=macro_bins[0][0], maxbinnum=macro_bins[-1][1])
    E_MIN, E_MAX = emin[0], emax[-1]

    logger.info('Checking PSF file...')
    PSF_FILE = kwargs['psffile']
    if PSF_FILE is not None:
        logger.info('Using %s'%PSF_FILE)
    else:
        logger.info('ATT: File not found: %s'%PSF_FILE)
        logger.info('Creating PSF file with gtpsf...')
        try:
            IRFS = kwargs['irfs']
        except:
            logger.info('ERROR: provide IRFS!')
            sys.exit()
        try:
            EVTYPE = kwargs['evtype']
        except:
            logger.info('ERROR: provide event type!')
            sys.exit()
        try:
            LT_FILE = data.LT_FILE
        except:
            LT_FILE = ''
        if os.path.exists(LT_FILE):
            logger.info('Livetime file found.')
        else:
            try:
                LT_FILE = kwargs['ltfile']
            except:
                logger.info('ERROR: provide livetime file or list!')
                sys.exit()
            if LT_FILE.lower().endswith(('.txt', '.dat')):
                lt_dict = {'infile1' : LT_FILE,
                           'outfile' : 'DEFAULT',
                           'chatter': 4,
                           'clobber': 'no'}
                from Xgam.utils.ScienceTools_ import gtltsum
                label_lt = IRFS + '_evt' + str(EVTYPE)
                out_gtltsum = gtltsum(label_lt,lt_dict)
                expcube = out_gtltsum
            else:
                expcube = LT_FILE
            
        from Xgam.utils.ScienceTools_ import gtpsf
        label_psf =  IRFS + '_evt' + str(EVTYPE)
        psf_dict = {'expcube' : expcube,
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
        PSF_FILE = gtpsf(label_psf, psf_dict)

    #NOTE: IMPLEMENT TYPE2 ??
    if kwargs['typesrcmask'] == 1:
        from Xgam.utils.mkmask_ import mask_src_fluxPSFweighted_1 as mask_src
    else:
        from Xgam.utils.mkmask_ import mask_src_fluxPSFweighted_2 as mask_src
    
    from Xgam.utils.mkmask_ import mask_gp
    from Xgam.utils.parsing_ import get_psf_en_univariatespline

    nside = kwargs['nside']
    out_label = kwargs['outflabel']
    npix = hp.nside2npix(nside)
    src_cat = kwargs['srccat']
    src_ext_cat = kwargs['srcextcat']
    out_name_list = os.path.join(X_OUT, 'fits/MaskSmart_'+out_label+'_list.txt')
    out_list = []

    for i, (minb, maxb) in enumerate(macro_bins):
        logger.info('Considering bins from %i to %i...' %(minb, maxb))
        emin, emax, emean = get_energy_from_fits(micro_bin_file, minbinnum=minb, maxbinnum=maxb)
        E_MIN, E_MAX = emin[0], emax[-1]
        E_MEAN = np.sqrt(emax[0]*emin[-1])
        logger.info('Energies %.1f - %.1f [MeV]' %(E_MIN, E_MAX))

        out_name = os.path.join(X_OUT, 'fits/MaskSmart_'+out_label+'_%.1f_MeV.fits'%(E_MIN))
        out_list.append(out_name)
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
            hp.write_map(out_name, mask, coord='G', overwrite=True)
            logger.info('Created %s \n' %out_name)

        if kwargs['show'] == True:
            import matplotlib.pyplot as plt
            hp.mollview(mask, cmap='bone')
            plt.show()

    logger.info('Writing list of output files: %s'%out_name_list)
    np.savetxt(out_name_list, out_list, fmt='%s')
    logger.info('Done!')

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkSmartMask(**args.__dict__)
