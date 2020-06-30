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


"""PolSpice functions
"""

import os
import sys
import numpy as np
import scipy as sp
from scipy.special import lpmv as leg_m_n
import healpy as hp
import pickle as pic
import astropy.io.fits as pf
from matplotlib import pyplot as plt

from Xgam import X_OUT
from Xgam.utils.logging_ import logger

def new_binning(xmin, xmax, nbin=25, bin_type='lin', out_type=int, custom_bins=None):
    """
    Define the new binning.

    Parameters
	----------

	Returns
	-------
	array
	   the array with the edges of the new binning
    """
    if bin_type == 'lin' and custom_bins is None:
        binning_ = np.linspace(xmin, xmax, num=nbin+1, dtype=out_type)
    elif bin_type == 'log' and custom_bins is None:
        if xmin == 0:
            xmin = 1
        binning_ = np.logspace(np.log10(xmin), np.log10(xmax), num=nbin+1, dtype=out_type)
    elif type(custom_bins) == list or type(custom_bins) == np.ndarray:
        binning_ = np.array(custom_bins)
    else:
        logger.info('ERROR: Invalid binning type. Choose lin or log, or customize it.')
        sys.exit()
        
    logger.info('Multipole binning:%s'%str(binning_))
    return binning_

def pol_ccf_parse(pol_ccf_out_file, pol_cov_out_file, rebin=None, pol='circular'):
    """
    Parser for the cross-correlation function.

    Parameters
	----------

	Returns
	-------
	array, array, array, tensor
    """
    hdu = pf.open(pol_cov_out_file)
    _cov = hdu[0].data[0]
    _cov = np.array(_cov)
    f = open(pol_ccf_out_file, 'r')

    th_, ccf_ = [], []
    for line in f:
        try:
            th, ccf = [float(item) for item in line.split()]
            th_.append(th)
            ccf_.append(ccf)
        except:
            pass
    th_ = np.array(th_)
    ccf_ = np.array(ccf_)

    #if rebin is not None:

    #Tranforming covariance
    #NOTE: think about Wl and Wpix correction!
    ccf_cov_ = np.zeros(np.shape(_cov))
    Pl_ = []
    if pol == 'circular':
        A=(2.*l+1.)/(4.*np.pi)
        leg_ord = 1
    elif  pol == 'dipolar':
        A=(2.*l+1.)/((l*(l+1.))*4.*np.pi)
        leg_ord = 2
    else:
        logging.info('ERROR: Invalid polarization type in covariance tranformation!')
        sys.exit()

    l_ = np.arange(len(_cov), dtype=int)
    for th in th_:
        Pl_.append(leg_m_n(leg_ord, l_, np.cos(np.radians(th))))
    Pl_ = np.array(Pl_)

    for i,th_1 in enumerate(th_):
        for j,th_2 in enumerate(th_):
            if j>=i:
                ccf_cov_[i,j] = np.sum((A[j]*Pl_[j])*np.transpose((A[i]*Pl_[i])*_cov))

    ccf_cov_ = ccf_cov_ + np.transpose(ccf_cov_)
    idx_diag = np.diag_indices(len(ccf_cov_))
    ccf_cov_[idx_diag] = ccf_cov_[idx_diag]/2.0

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
    pol_config = config_file_name
    pol_config_file = open(pol_config, 'w')
    for key in pol_dict:
        pol_config_file.write('%s = %s \n'%(key, str(pol_dict[key])))
    logger.info(config_file_name)
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

def pol_cl_parse(pol_cl_out_file, pol_cov_out_file, wl_array=None, rebin=None, nbin=25, 
                        lmin=1, lmax=1500, bin_type='lin', custom_bins=None, show=False):
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
    if wl_array is not None:
        wl = wl_array
        _l = _l[:len(wl)]
        _cl = (_cl[:len(wl)])/(wl**2)
        _cov_ = np.array([_cov[i][:len(wl)] for i in range(0,len(wl))])
        _cov_ = _cov_/(wl**2)
        for l in range(len(wl)):
            _cov_[l] = _cov_[l]/(wl**2)
    else:
        pass

    if rebin:
        logger.info('Rebinning in multipole range')
        _lr, _clr, _clerrr = [], [], []
        if custom_bins is not None:
            logger.info('Custom multipole binning.')
            rebinning = new_binning(lmin, lmax, nbin, bin_type=bin_type, custom_bins=custom_bins)  
        else:
            rebinning = new_binning(lmin, lmax, nbin, bin_type=bin_type)
        for bmin, bmax in zip(rebinning[:-1], rebinning[1:]):
            logger.info('considering %i < li < %i'%(bmin,bmax))
            _index = np.where(np.logical_and(_l>=bmin, _l<bmax))[0]
            _index = _index.astype(int)
            if bin_type=='log':
                _lmean = np.sqrt(bmin*bmax)
            if bin_type=='lin':
                _lmean = (bmin+bmax)/2
            _lr.append(_lmean)
            _clmean = np.mean(_cl[_index])
            _clerr = np.mean(_cov[bmin:bmax,bmin:bmax])
            logger.info('cl_mean %.3e'%_clmean)
            logger.info('cl_mean err %.3e'%np.sqrt(_clerr))
            _clr.append(_clmean)
            _clerrr.append(np.sqrt(_clerr))
        _l = np.array(_lr)
        _cl = np.array(_clr)
        _clerr = np.array(_clerrr)
    else:
        _clerr = np.array([np.sqrt(_cov[i][i]) for i in _l.astype(int)])
        
    if show:
        plt.figure()
        plt.errorbar(_l, _cl, yerr=_clerr, fmt='.')
        plt.plot([_l[0], _l[-1]], [0, 0], c='silver', linewidth=1)
        plt.title('CAPS', size=20)
        plt.xlabel('$\ell$', size=18)
        plt.ylabel('C$_{\ell}$', size=18)
        plt.xscale('log')
        plt.show()
     
    return np.array(_l), np.array(_cl),  np.array(_clerr)

def pol_cov_parse(pol_cov_out_file, wl_array=None, rebin=None, nbin=25, bin_type='lin', 
                  lmin=1, lmax=1500, custom_bins=None, show=False):
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
    hdu = pf.open(pol_cov_out_file)
    _cov = hdu[0].data[0]
    hdu.close()
    _l = np.arange(len(_cov))
    if wl_array is not None:
        wl = wl_array
        _l = np.arange(len(wl_array))
        _cov = np.array([_cov[i][:len(wl)] for i in range(0,len(wl))])
        _cov = _cov/(wl**2)
        for l in _l:
            _cov[l] = _cov[l]/(wl[l]**2)
    if rebin:
        _covr = []
        _lr = []
        if custom_bins is None:
            rebinning = new_binning(lmin, lmax, nbin, bin_type=bin_type)
        else:
            rebinning = new_binning(lmin, lmax, nbin, bin_type=bin_type, custom_bins=custom_bins)  
        for imin, imax in zip(rebinning[:-1], rebinning[1:]):
            _imean = np.sqrt(imin*imax)
            _covrj = []
            if bin_type=='log':
                _lmean = np.sqrt(imin*imax)
            if bin_type=='lin':
                _lmean = (imin+imax)/2
            _lr.append(_lmean)
            for jmin, jmax in zip(rebinning[:-1], rebinning[1:]):
                 _covrj.append(np.mean(_cov[imin:imax, jmin:jmax]))
            _covr.append(np.array(_covrj))
        _cov = np.array(_covr)
        _l = np.array(_lr)
    else:
        pass
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
        cax = ax.matshow(_cov2ploti, origin='lower', 
                        extent=[_l[0], _l[-1], _l[0], _l[-1]], 
                        aspect='auto', cmap='viridis')
        plt.title('$\sigma_{ij}/\sqrt{\sigma_{ii}\sigma_{jj}}$')
        plt.xlabel('$l_{i}$')
        plt.ylabel('$l_{j}$')
        cb = plt.colorbar(cax, format='$%i$')
        plt.grid()
        plt.show()
        
    return _cov

def pol_cl_calculation(pol_dict, config_file_name, wl_array=None, rebin=None, lmin=1, 
                      lmax=1500, nbin=25, custom_bins=None, bin_type='lin', show=False):
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
    wl = wl_array
    r = rebin
    s = show
    nb = nbin
    bt = bin_type
    cbins = custom_bins
    lm, lM = lmin, lmax
    pol_cl_out_file = pol_dict['clfile']
    pol_cov_out_file = pol_dict['covfileout']
    config_file = pol_create_config(pol_dict, config_file_name)
    pol_run(config_file)
    
    if custom_bins is None:
        _l, _cl, _clerr = pol_cl_parse(pol_cl_out_file, pol_cov_out_file, wl_array=wl, 
                                         rebin=r, nbin=nb, lmin=lm, lmax=lM, bin_type=bt)
        _cov = pol_cov_parse(pol_cov_out_file,  wl_array=wl, lmin=lm, lmax=lM, rebin=r,  
                                                            nbin=nb, bin_type=bt, show=s)
    else:
        _l, _cl, _clerr = pol_cl_parse(pol_cl_out_file, pol_cov_out_file, lmin=lm, lmax=lM,
                            wl_array=wl, rebin=r, nbin=nb, bin_type=bt, custom_bins=cbins)  
        _cov = pol_cov_parse(pol_cov_out_file,  wl_array=wl, rebin=r, lmin=lm, lmax=lM,
                                        nbin=nb, bin_type=bt,  custom_bins=cbins, show=s)                           
    
    return _l, _cl, _clerr, _cov

def main():
    """test module
    """
    return 0


if __name__ == '__main__':
    main()
