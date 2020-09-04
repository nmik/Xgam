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

"""Simple model fit of the APS
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
PARSER.add_argument('--outflabel', type=str, required=True,
                    help='Output file label')
PARSER.add_argument('--lminfit', type=int, default=10,
                    help='minimum multippole to consider in the fit')
PARSER.add_argument('--lmaxfit', type=int, default=1000,
                    help='maximum multippole to consider in the fit')
PARSER.add_argument('--title', type=str, default='Cross-correlation', 
                    help='Plot title')
PARSER.add_argument('--xscale', type=str, choices=['log', 'linear'], default='linear', 
                    help='x-axis scale')
PARSER.add_argument('-m', '--model', type=str, choices=['1h', '1h2h'], default='1h', 
                    help='Fit model to use')
    
def onehalo_model(x, c):

    return c
    
def twohalo_model(x, a, b, c):

    return a * x**(-b) + c
    
    
def aps_fit(**kwargs):
    """
    Viewer interface for APS and Covariance matrix                                                                                                                                                             
    """
    
    params = parse_polspice_aps(kwargs['infile'])
    _emin, _emax = params[0], params[1]
    _l, _cl, _clerr = params[2], params[3], params[4]
    _cov_ = params[5]
	
    logger.info('Generating APS plots...')
    for i, cl in enumerate(_cl):
        valid = ~(np.isnan(cl))
        l = _l[i][valid]
        cl = cl[valid]
        clerr = _clerr[i][valid]
        fitrange = np.where((l>int(kwargs['lminfit']))&(l<int(kwargs['lmaxfit'])))[0]
	    
        fig = plt.figure()
        plt.plot([0, 2000], [0, 0], '--', color='silver')
        aps_label = '%.1f-%.1f GeV'%(_emin[i]/1000, _emax[i]/1000)
        plt.errorbar(l, cl, fmt='o', markersize=5, elinewidth=2, label='Data', 
                     color='0.3', yerr=clerr)
        if kwargs['model'] == '1h':
            popt, pcov = curve_fit(onehalo_model, l[fitrange], cl[fitrange], 
                                         sigma=clerr[fitrange], absolute_sigma=True) 
            print('-----------------------')
            print('Fit range: %i - %i'%(kwargs['lminfit'], kwargs['lmaxfit']))
            print('Best-fit Cp = %.2e+-%2e'%(popt[0], np.sqrt(pcov[0][0])))
            print('Detection Significance: %.3f'%(popt[0]/np.sqrt(pcov[0][0])))
            print('-----------------------')
            
            plt.fill_between(l[fitrange], np.full(len(l[fitrange]), popt[0]-np.sqrt(pcov[0][0])), 
                             np.full(len(l[fitrange]), popt[0]+ np.sqrt(pcov[0][0])), 
                             color='lightcoral', alpha=0.5)
            plt.plot([l[fitrange][0], kwargs['lmaxfit']], [popt[0], popt[0]], '--', color='red', label='C$_{p}$ = %.1e'%popt[0])
        else:
            popt, pcov = curve_fit(twohalo_model, l[fitrange], cl[fitrange], 
                                         sigma=clerr[fitrange], absolute_sigma=True) 
                                         
            chi2 = np.sum(((cl[fitrange] - twohalo_model(l[fitrange], *popt))/clerr[fitrange])**2)
            chi2red = chi2/len(fitrange)
            
            print('-----------------------')
            print('Fit range: %i - %i'%(kwargs['lminfit'], kwargs['lmaxfit']))
            print('Best-fit params = %s'%(popt))
            print('Chi2red: %.3f'%(chi2red))
            print('-----------------------')
            
            plt.plot(l[fitrange], twohalo_model(l[fitrange], *popt), '--', 
                     color='red', label='Best fit: %.1e $\ell^{-%.1f}$ + (%.1e)'%(popt[0], popt[1], popt[2]))
	        
        plt.xlim(kwargs['lminfit']-10, kwargs['lmaxfit']+5)
        #plt.ylim(-1e-11,1e-11)
        plt.xscale(kwargs['xscale'])
        plt.xlabel('Multipole', size=15)
        plt.ylabel('C$^{sig}_{\ell}$', size=15) #[(cm$^{-2}$s$^{-1}$sr$^{-1}$)sr]
        plt.title(kwargs['title'])
        plt.tight_layout()
        plt.legend(loc=1, fontsize=16)
        
        plt.savefig('output/figs/%s_fit_%s.png'%(kwargs['outflabel'], kwargs['model']))
    
    plt.show()


	
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    aps_fit(**args.__dict__)