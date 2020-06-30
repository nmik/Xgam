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

"""Visualize the healpix maps
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

from Xgam import X_OUT
from Xgam.utils.matplotlib_ import *
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.parsing_ import parse_datafluxmaps

__description__ = 'Visualize the maps'


"""Command-line switches.
"""
formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infile', type=str, required=True,
                    help='input fits file')
PARSER.add_argument('--title', type=str, default=' ', 
                    help='map title')


E_IGRB_sub = [123.57,176.22,247.82,346.12,473.48,699.28,983.44,1344.94,1944.00,2696.27,
              3739.12,5483.63,7767.11,11083.96,15803.47,22854.50,31694.83,44270.961,
              64014.963,88160.632,125660.80,175498.30,246721.16,351871.07,473879.41,679459.47]
E2F_IGRB_sub = [0.0006370399,0.000530884,0.00045098,0.000410,0.000394,0.000394,0.000341,
                0.000279,0.0001795,0.000161,0.000132,0.000125,0.000122,0.000166,
                0.000127,0.000115,0.0000962,0.0000917,0.0000764,0.0000637,0.0000398,
                0.0000348,.0000237,0.0000218,0.00000688,0.00000215]		  
E_IGRB_up	= [121.94,175.12,248.02,351.26,480.54,709.72,977.50,1374.73,2072.83,2681.72,
               4043.40,5530.29,7943.05,11017.75,16495.20,22561.64,32397.78,61454.28,61454.28,
               89478.77,130230.12,179351.34,245226.01,344841.64,505033.64,699716.11]
E2F_IGRB_up = [0.0009623,0.00084950,0.00077923, 0.0007011,.000707945,0.00070794,
                0.00062493,0.0005463,0.0004823,0.00045098,0.00039054,0.00034145,
                0.00033176,0.0003162,0.0002560,0.0002282,0.0001920,0.00014399,0.00014399,
                0.00011659,0.0000721,0.0000601,0.0000425,0.0000348,0.00001333,0.00000589]

def ref_igrb_band():
    """
    Draw the systematics band of the Pass7 IGRB flux
    """

    up_line_y2 = np.interp(E_IGRB_sub, E_IGRB_up, E2F_IGRB_up)
    lab='Ackermann et al. 2015'
    igrb = plt.fill_between(E_IGRB_sub, E2F_IGRB_sub, up_line_y2, label=lab,
                    facecolor='none', hatch='//', edgecolor='black', zorder=10)
    return igrb, lab


def flux_view(**kwargs):
    """
    Viewer interface for flux spectrum                                                                                                                                                             
    """
    
    params = parse_datafluxmaps(kwargs['infile'])
    _emin, _emax, _emean = params[0], params[1], params[2]

    if len(params[-1]) == 0:
        _f, _ferr = params[3], params[4]

    else:
        _n, _nsx, _ndx = params[-6], params[-5], params[-4]
        _f, _fsxerr, _fdxerr = params[-3], params[-2], params[-1]
        _ferr = [abs(_f-_fsxerr), abs(_f-_fdxerr)]
	
    logger.info('Plotting Energy Spectrum ...')
    fig, ax = plt.subplots(1)
    leg, lab = ref_igrb_band()

    spec = plt.errorbar(_emean, (_f)/(_emax-_emin)*_emean**2, fmt='o', color='red',
                        markersize=3, elinewidth=2, xerr=[(_emean-_emin), (_emax-_emean)], 
                        yerr=_ferr/(_emax-_emin)*_emean**2)
                        
    Aeff_sys = plt.fill_between(_emean, (_f-0.1*_f)/(_emax-_emin)*_emean**2, 
                                (_f+0.1*_f)/(_emax-_emin)*_emean**2, color='red', alpha=0.3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('Energy [MeV]', size=14)
    plt.ylabel('E$^{2}$ $\cdot$ Flux [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]', size=14)
    plt.title(kwargs['title'])
    plt.legend([leg, spec, Aeff_sys], [lab, 'UGRB flux (Xgam)', 'A$_{eff}$ sys. band'], loc=3, fontsize=15)

    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim(1e-7, 4.e-2)
    plt.xlim(100, 1000000)
    
    if len(params[-1]) != 0:
        logger.info('Plotting Foreground best-fit Normalization ...')
        
        plt.figure()
        plt.title(kwargs['title'])
        plt.fill_between(_emean, _nsx, _ndx, color='lightcoral', alpha=0.5)
        plt.plot(_emean, _n, 'ro-')
        plt.xscale('log')
        plt.xlabel('Energy [MeV]', size=14)
        plt.ylabel('IEM Normalization', size=14)
	    
    plt.show()

	
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    flux_view(**args.__dict__)