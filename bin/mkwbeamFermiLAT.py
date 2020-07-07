#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


"""Analysis module                                                              
"""


import os
import ast
import argparse
import numpy as np
import healpy as hp
from itertools import combinations

__description__ = 'Makes the Fermi LAT Wbeam functions'


"""Command-line switches.                                                       
"""


from Xgam import X_OUT, FT_DATA_FOLDER
from Xgam.utils.wbeamfunc_ import * 
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.parsing_ import get_psf_th_en_bivariatespline

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--psffile', type=str, required=False,
        help='PSF input file (if given --irfs and --evtype are discarded)', 
        default=None)
PARSER.add_argument('--irfs', type=str, required=False,
                    help='The Fermi-LAT IRFs to consider')
PARSER.add_argument('--evtype', type=int, required=False,
                    help='The event type to consider')
PARSER.add_argument('--ltfile', type=str, required=False,
                    help='livetimes .fits or .txt with list of .fits', 
                    default=None)
PARSER.add_argument('--lmax', type=str, required=False,
                    help='Maximuma multipole value', 
                    default=1500)
PARSER.add_argument('--show', type=ast.literal_eval, choices=[True, False],
                    help='True if you want to visualize the Wbeam matrix',
                    default=False)

def mkwbeamFermiLAT(**kwargs):
    """ 
    This routine generates the W_beam functions from Fermi-LAT 
    Instrument Response Functions (IRFs)
                                       
    """
    logger.info('Checking PSF file...')
    PSF_FILE = kwargs['psffile']
    if PSF_FILE is not None:
        logger.info('Using %s'%PSF_FILE)
    else:
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
                out_gtltsum = gtltsum(label_lt, lt_dict)
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
                    'emin' : 100,
                    'emax' : 1000000,
                    'nenergies' : 500,
                    'thetamax' : 30.0,
                    'ntheta' : 100,
                    'chatter': 4,
                    'clobber': 'no'}
        PSF_FILE = gtpsf(label_psf, psf_dict)
              
    logger.info('Calculating Wbeam Function...')
    if PSF_FILE is not None:
    	out_wb_txt = os.path.join(X_OUT, os.path.basename(PSF_FILE).replace('_psf.fits', 
    	                                  '_Wbeam.txt'))
    else:
    	out_wb_txt = os.path.join(X_OUT, '%s_evt%i_Wbeam.txt'%(IRFS, EVTYPE))

    psf = get_psf_th_en_bivariatespline(PSF_FILE)
    ELL_MAX = kwargs['lmax']
    l_ = np.arange(ELL_MAX)
    wb = build_wbeam(psf, l_, out_wb_txt)
    
    logger.info('Created %s'%out_wb_txt)

    if kwargs['show']==True:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(wb, origin='lower', aspect='auto', cmap='viridis')
        ax.xaxis.set_ticks_position('bottom')
        plt.xlabel('Multipole $l$', size=15)
        plt.ylabel('Energy', size=15)
        cbar = plt.colorbar(cax)
        cbar.set_label('$W_{beam}(l, E)$', size=15)
        plt.grid(False)
        
        plt.show() 

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkwbeamFermiLAT(**args.__dict__)
