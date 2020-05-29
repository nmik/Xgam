#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


"""PolSpice functions
"""

import os
import numpy as np
import scipy as sp
import healpy as hp
import pickle as pic
import astropy.io.fits as pf
from matplotlib import pyplot as plt

from Xgam import X_OUT
from Xgam.utils.logging_ import logger

def new_binning():
    """
    Define the new binning.
    
    Parameters
	----------
	
	Returns
	-------
	array
	   the array with the edges of the new binning
    """
    binning_ = np.array([]])
    return binning_
    
def pol_ccf_parse(pol_ccf_out_file, pol_cov_out_file, rebin=None):
    """
    Parser for the cross-correlation function.
    
    Parameters
	----------
	
	Returns
	-------
	array, array, array, tensor
    """
    return th_, ccf_, ccferr_, ccf_cov_
    
def pol_create_config(pol_dict, config_file_name):
    """
    Creates and returns PolSpice config ascii file.
    
    Parameters
	----------   
    pol_dict : python dict
        a dictionary where all the parameters of a tipical PolSpice config
        ascii file should have.
    config_file_name : str
        name of the confg ascii file that will be created
          
    Returns
    -------
    str
        Name of the created PolSpice config file
    """
    if not os.path.exists(X_OUT+'output_polspice'):
    os.makedirs(X_OUTX_OUT+'output_polspice')
    pol_config = os.path.join(X_OUT, 'output_polspice/%s.txt'%config_file_name)
    pol_config_file = open(pol_config, 'w')
    for key in pol_dict:
        pol_config_file.write('%s = %s \n'%(key, str(pol_dict[key])))
    logger.info('Created output/output_polspice/%s.txt'%config_file_name)
    return pol_config 

def pol_run(config_file):
    """
    Runs PolSpice.
    
    Parameters
	----------  
	config_file : str 
          configuration txt file
    """
    os.system('spice -optinfile %s' %(config_file))
    
    return 0

def pol_cl_parse(pol_cl_out_file, pol_cov_out_file, raw_corr=None, rebin=None):
    """
    Parser of the txt output file of PolSpice, which contains the angular power spectrum (APS).
    
    Parameters
	----------
    pol_cl_out_file : str
        ascii output file created by PolSpice 
    pol_cov_out_file : str
        .fits file containing the covariance matrix created by PolSpice
    raw_corr : (float, numpy array)
        must be a python list with 2 entries: the first one must be the
        white poissonian noise, the second must be the array (or the spline)
        of the  Wbeam function as a funcion of l integrated in a energy bin.
    rebin : list (or array)
        if not None, it has to be the list defining the edges of the new binning.
    
    Returns
    -------
    array, array, array
        in order: array of the ells, array of the Cls, array of the Cls errors 
        (estimated from the covariance matrix). If rebin is not None the arrays are 
        the rebinned array. 
    """
    hdu = pf.open(pol_cov_out_file)
    _cov = hdu[0].data[0]
    _invcov = np.linalg.inv(_cov)    
    f = open(pol_cl_out_file, 'r')

    _l, _cl = [], []
    for line in f:
        try:
            l, cl = [float(item) for item in line.split()]
            _l.append(l)
            _cl.append(cl)
        except:
            pass
    _l = np.array(_l)
    _cl = np.array(_cl)
    if raw_corr is not None:
        cn, wl = raw_corr[0], raw_corr[1]
        _l = _l[:len(wl)]
        _cl = (_cl[:len(wl)] - cn)/(wl**2)
        _cov = np.array([_cov[i][:len(wl)] for i in range(0,len(wl))])
        _cov = _cov/(wl**2)
        for l in _l:
            _cov[l] = _cov[l]/(wl[l]**2)
    else:
        pass
# TO BE REWRITTEN --------------------------------------------------
#    if rebin:
#         logger.info('Rebinning in multipole range:')
#         logger.info('%s'%str(rebinning))
#         _lr, _clr, _clerrr = [], [], []
#         for bmin, bmax in zip(rebinning[:-1], rebinning[1:]):
#             logger.info('considering %i < li < %i'%(bmin,bmax))
#             _index = np.where(np.logical_and(_l>=bmin, _l<bmax))[0]
#             _index = _index.astype(int)
#             _lmean = np.sqrt(bmin*bmax)
#             _lr.append(_lmean)
#             _clmean = np.mean(_cl[_index])
#             _clerr = np.mean(_cov[bmin:bmax,bmin:bmax])
#             logger.info('cl_mean %.3f'%_clmean)
#             logger.info('cl_mean err %.3f'%np.sqrt(_clerr))
#             _clr.append(_clmean)
#             _clerrr.append(np.sqrt(_clerr))
#         _l = np.array(_lr) 
#         _cl = np.array(_clr)
#         _clerr = np.array(_clerrr)
#     else:
#         _clerr = np.array([np.sqrt(_cov[i][i]) for i in _l.astype(int)])
    return np.array(_l), np.array(_cl),  np.array(_clerr)

def pol_cov_parse(pol_cov_out_file, wl_array=None, rebin=None, show=False):
    """
    Parser of the fits output file of PolSpice containing the covariance 
    matrix of the angular power spectra (APS).

    Parameters
	----------
    pol_cov_out_file :
        .fits file containing the covariance matrix created by PolSpice
    wl_array : numpy array (or spline)
        array (or the spline) of the  Wbeam function as a funcion of l 
        integrated in a energy bin.
    rebin : list (or array)
        if not None, it has to be the list defining the edges of the new binning. 
    show : bool
        if True the plot of the covariance matrix is shown
          
    Returns
    -------
    numpy tensor
        rebinned covariance matix (shaped as a tensor NxN, where N is the number of bins)    
    """
    show = False
    hdu = pf.open(pol_cov_out_file)
    _cov = hdu[0].data[0]
    hdu.close()
    _l = np.arange(len(_cov))
    if wl_array is not None:
        wl = wl_array
        _l = np.arange(len( wl_array))
        _cov = np.array([_cov[i][:len(wl)] for i in range(0,len(wl))])
        _cov = _cov/(wl**2)
        for l in _l:
            _cov[l] = _cov[l]/(wl[l]**2)   
# TO BE REWRITTEN -------------------------------------------------- 
#     if rebin:
#         _covr = []
#         _lr = []
#         for imin, imax in zip(rebinning[:-1], rebinning[1:]):
#             _imean = np.sqrt(imin*imax)
#             _covrj = []
#             _lmean = np.sqrt(imin*imax)
#             _lr.append(_lmean)
#             for jmin, jmax in zip(rebinning[:-1], rebinning[1:]):
#                 _covrj.append(np.mean(_cov[imin:imax, jmin:jmax]))
#             _covr.append(np.array(_covrj))
#         _cov = np.array(_covr)
#         _l = np.array(_lr)
#     else:
#         pass 
    pic.dump(_cov, 
             open(pol_cov_out_file.replace('.fits', '.pkl'),'wb'))
    if show==True:
         _cov2ploti = []
         for i in range(0, len(_l)):
             sigii = _cov[i][i]
             _cov2plotj = []
             for j in range(0, len(_l)):
                 sigjj = _cov[j][j]
                 sigij = _cov[j][i]
                 if sigij < 0:
                     sigij = 1e-100
                 _cov2plotj.append(np.sqrt(sigij/np.sqrt(sigii*sigjj)))
             _cov2ploti.append(_cov2plotj)
         _cov2ploti = np.array(_cov2ploti)
         fig = plt.figure(facecolor='white')
         ax = fig.add_subplot(111)
         cax = ax.matshow(np.log10(np.abs(_cov)), origin='lower', 
                          aspect='auto', cmap='Spectral')
         en_tick = list(np.logspace(0, np.log10(1500), 6).astype(int))
         ax.set_yticklabels(['']+en_tick)
         ax.set_xticklabels(['']+en_tick)
         plt.title('Covariance matrix')
         plt.xlabel('$l_{i}$')
         plt.ylabel('$l_{j}$')
         cb = plt.colorbar(cax, format='$%i$')
         plt.grid()
         save_current_figure(os.path.basename(pol_cov_out_file).replace('.fits', ''))
    return _cov

def pol_cl_calculation(pol_dict, config_file_name, raw_corr=None, rebin=None, show=False):
    """ 
    Creates and runs the PolSpice config file, and returns the arrays of
    1) the multipoles (rebinned or not); 2) the Cl corresponding to the 
    multipoles; 3) the errors associated to the Cls; 4) the covariance 
    matrix.

    Parameters
	----------
    pol_dict : python dict
        a dictionary where all the parameters of a tipical PolSpice config
        ascii file should have.
    config_file_name : str
        name of the confg ascii file that will be created
    raw_corr : [float, numpy array]
        must be a python list with 2 entries: the first one must be the
        white poissonian noise, the second must be the array (or the spline)
        of the  Wbeam function as a funcion of l integrated in a energy bin.
    rebin : None
        if True a multipole rebinni of the APS is done.
    show : bool
        if True the plot of the covariance matrix is shown
        
    Returns
    -------
    array, array, array, tensor
        in order: array of the ells, array of the Cls, array of the error on the Cls
        errors, tensor of the covariance matrix. If rebin is not None the arrays are 
        the rebinned array. 
    """
    corr = raw_corr
    r = rebin
    s = show
    pol_cl_out_file = pol_dict['clfile']
    pol_cov_out_file = pol_dict['covfileout']
    config_file = pol_create_config(pol_dict, config_file_name)
    pol_run(config_file)
    cn, wl = 0., None
    if raw_corr is not None:
        cn, wl = raw_corr[0], raw_corr[1]
    _l, _cl, _clerr= pol_cl_parse(pol_cl_out_file, pol_cov_out_file,
                                  raw_corr=corr,rebin=r)
    _cov = pol_cov_parse(pol_cov_out_file,  wl_array=wl, 
                         rebin=r, show=s)
    return _l, _cl, _clerr, _cov
    
def main():
    """test module
    """
    print rebinning
    


if __name__ == '__main__':
    main()
