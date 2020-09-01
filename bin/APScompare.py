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

"""Compare APS 
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
PARSER.add_argument('-f', '--infile', type=str, nargs='*', required=True,
                    help='input txt files output of mkcrosscorrelation.py')
PARSER.add_argument('-l', '--labels', type=str, nargs='*', default=None,
                    help='input txt files output of mkcrosscorrelation.py')
PARSER.add_argument('--title', type=str, default='Cross-correlation', 
                    help='Plot title')
PARSER.add_argument('--xscale', type=str, choices=['log', 'linear'], default='linear', 
                    help='x-axis scale')
    
def aps_compare(**kwargs):
    """
    Viewer interface for APS and Covariance matrix                                                                                                                                                             
    """
    
    logger.info('Generating comparison APS plots...')
    
    l_f, cl_f, clerr_f = [], [], []
    
    for f in kwargs['infile']:
        params = parse_polspice_aps(f)
        _emin, _emax = params[0], params[1]
        _l, _cl, _clerr = params[2], params[3], params[4]
        _cov_ = params[5]
        
        l_f.append(_l)
        cl_f.append(_cl)
        clerr_f.append(_clerr)
        
    l_f = np.array(l_f)
    cl_f = np.array(cl_f)
    clerr_f = np.array(clerr_f)
    
    nfiles = l_f.shape[0]
    nebins = l_f.shape[1]
    
    lab = kwargs['labels']
    if kwargs['labels'] is None:
        lab = []
        for i in range(nfiles):
            lab.append('File %s'%i)
    
    for i in range(nebins):
        plt.figure()
        
        plt.hlines(0, 1, 1000, color='silver')
        
        plt.title('%s'%kwargs['title'])
        for j in range(nfiles):
            plt.errorbar(l_f[j][i], cl_f[j][i], fmt='o', markersize=5, elinewidth=2,
                         yerr=clerr_f[j][i], label=lab[j])
        
        plt.xscale(kwargs['xscale'])
        plt.xlabel('$\ell$', size=15)
        plt.ylabel('C$_{\ell}$', size=15)
        
        plt.legend()
   
    plt.show()


	
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    aps_compare(**args.__dict__)