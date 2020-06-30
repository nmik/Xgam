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
from Xgam.utils.parsing_ import parse_polspice_aps

__description__ = 'Visualize the maps'


"""Command-line switches.
"""
formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infile', type=str, required=True,
                    help='input fits file')
PARSER.add_argument('--title', type=str, default='Cross-correlation', 
                    help='Plot title')
PARSER.add_argument('--xscale', type=str, choices=['log', 'linear'], default='linear', 
                    help='x-axis scale')

def aps_view(**kwargs):
    """
    Viewer interface for APS and Covariance matrix                                                                                                                                                             
    """
    
    params = parse_polspice_aps(kwargs['infile'])
    _emin, _emax = params[0], params[1]
    _l, _cl, _clerr = params[2], params[3], params[4]
    _cov_ = params[5]
	
    logger.info('Generating APS plots...')
    for i, cl in enumerate(_cl):
    
        sig = cl/_clerr[i]
        print('-----------------------')
        print('Detection Significance:')
        print(sig)
        print('-----------------------')
	    
        aps_label = '%.1f-%.1f GeV'%(_emin[i]/1000, _emax[i]/1000)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot([0, 2000], [0, 0], '--', color='silver')
        plt.errorbar(_l[i], cl, fmt='o', markersize=5, elinewidth=2, label=aps_label, 
                     color='red', yerr=_clerr[i])
        plt.xlim(np.amin(_l[i])-1, np.amax(_l[i])+1)
        plt.xscale(kwargs['xscale'])
        plt.xlabel('Multipole', size=15)
        plt.ylabel('C$^{sig}_{\ell}$ [(cm$^{-2}$s$^{-1}$sr$^{-1}$)sr]', size=15)
        plt.title(kwargs['title'])
        plt.tight_layout()
        plt.legend(loc=1, fontsize=16)
        
    logger.info('Generating Covariance matrices plots...')
    for k, cov in enumerate(_cov_):
    
        _cov2plotij_ = []
        for i in range(0, len(cov)-1):
            sigii = cov[i][i]
            _cov2plotj = []
            for j in range(0, len(cov[0])-1):
                sigjj = cov[j][j]
                sigij = cov[j][i]
                if sigij < 0:
                    sigij = 1e-100
                _cov2plotj.append(np.sqrt(sigij/np.sqrt(sigii*sigjj)))
            _cov2plotij_.append(_cov2plotj)
        _cov2plotij_ = np.array(_cov2plotij_)
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(_cov2plotij_, origin='lower', aspect='auto', cmap='viridis', 
                         extent=[_l[k][0], _l[k][-1], _l[k][0], _l[k][-1]])
        ax.xaxis.set_ticks_position('bottom')
        plt.title(kwargs['title']+' (%.1f-%.1f GeV)'%(_emin[k]/1000, _emax[k]/1000))
        plt.xlabel('Multipole $\ell_{i}$', size=15)
        plt.ylabel('Multipole $\ell_{j}$', size=15)
        cbar = plt.colorbar(cax)
        cbar.set_label('$\sigma_{ij}/\sqrt{\sigma_{ii}\sigma_{jj}}$', size=15)
        #plt.tight_layout()
        plt.grid(False)
        
    
    plt.show()


	
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    aps_view(**args.__dict__)