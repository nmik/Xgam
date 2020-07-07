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

import os
import re
import numpy as np
import healpy as hp
from numba import jit
import astropy.io.fits as pf
from itertools import product
from scipy.special import factorial
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.spline_ import xInterpolatedUnivariateSplineLinear

FORE_EN = re.compile('\_\d+\.')

def get_fore_integral_flux_map(fore_files_list, e_min, e_max):
    """
    
    A powerlaw is assumed for the foreground energy spectrum, hence
    the interpolation between 2 given maps at given energies (given
    by the model) is done in logarithmic scales.

    Parameters
    ----------
    fore_files_list: list of str
           Ordered list of the foreground files (one for each energy)
    e_min: float
           the min of the energy bin
    e_max: float
           the max of the energy bin
           
    Returns
    -------
    array
        foreground map integrated between e_min and e_max
        
    """
    fore_en = []
    for ff in fore_files_list:
        m = re.search(FORE_EN, ff)
        en = int(m.group(0).replace('_', '').replace('.', ''))
        fore_en.append(en)
    fore_en = np.array(fore_en)
    out_name = fore_files_list[0].replace('_%i.fits'%fore_en[0],
                                          '_%d-%d.fits'%(e_min, e_max))
    if os.path.exists(out_name):
        logger.info('ATT: file %s already exists and returned...'%out_name)
        fore_map = hp.read_map(out_name)
        return fore_map
    else:
        logger.info('Computing the integral flux of the foreground model...')
        logger.info('...between %.2f - %.2f'%(e_min, e_max))
        fore_emin_sx, fore_emin_dx = find_outer_energies(e_min, fore_en)
        fore_emax_sx, fore_emax_dx = find_outer_energies(e_max, fore_en)
        fore_emin_sx_ind = np.where(fore_en == fore_emin_sx)[0][0]
        fore_emin_dx_ind = np.where(fore_en == fore_emin_dx)[0][0]
        fore_emax_sx_ind = np.where(fore_en == fore_emax_sx)[0][0]
        fore_emax_dx_ind = np.where(fore_en == fore_emax_dx)[0][0]
        fore_fmin_sx = hp.read_map(fore_files_list[fore_emin_sx_ind])
        fore_fmin_dx = hp.read_map(fore_files_list[fore_emin_dx_ind])
        fore_fmax_sx = hp.read_map(fore_files_list[fore_emax_sx_ind])
        fore_fmax_dx = hp.read_map(fore_files_list[fore_emax_dx_ind])
        m1 = (np.log10(fore_fmin_sx)-np.log10(fore_fmin_dx))/ \
            (np.log10(fore_emin_sx)-np.log10(fore_emin_dx))
        m2 = (np.log10(fore_fmax_sx)-np.log10(fore_fmax_dx))/ \
            (np.log10(fore_emax_sx)-np.log10(fore_emax_dx))
        logfore1 = m1*(np.log10(e_min)-np.log10(fore_emin_sx))+ \
            np.log10(fore_fmin_sx)
        logfore2 = m2*(np.log10(e_max)-np.log10(fore_emax_sx))+ \
            np.log10(fore_fmax_sx)
        fore1 = 10**(logfore1)
        fore2 = 10**(logfore2)
        fore_integ_map = np.sqrt(fore1*fore2)*(e_max - e_min)
        hp.write_map(out_name, fore_integ_map)
        logger.info('Created file %s'%out_name)
        return fore_integ_map

def find_outer_energies(en_val, en_arr):
    """
    Gives the first element on the right and the first on the left
    of a given E value (en_val), among all values in an ordered array (en_arr).
    These values are used to integrate the foreground model in the considered energy bin.

    Parameters
    ----------
    en_val : float
           mean energy
    en_arr : float
           array of the energies at which the foreground model is given.
           
    Returns
    -------
    float, float
        first element on the right and the first on the left
    
    """
    en_sx_arr = en_arr[en_arr < en_val]
    en_dx_arr = en_arr[en_arr > en_val]
    if en_sx_arr.size == 0:
        en_sx = en_dx_arr[0]
        en_dx = en_dx_arr[1]
        logger.info('Considering model in the interval %.1f-%.1f MeV'
                    %(en_sx, en_dx))
    elif en_dx_arr.size == 0:
        en_sx = en_sx_arr[-2]
        en_dx = en_sx_arr[-1]
        logger.info('Considering model in the interval %.1f-%.1f MeV'
                    %(en_sx, en_dx))
    else:
        en_sx = en_sx_arr[-1]
        en_dx = en_dx_arr[0]
    return en_sx, en_dx

@jit
def myfactorial(_x):
    _facx = []
    for x in _x:
        n = 1
        for i in range(2, x+1):
            n *= i
        _facx.append(n)
    return np.array(_facx)
    
@jit
def poisson_likelihood(norm_guess, const_guess, fore_map, data_map, exp=None, sr=None):
    """
    Compute the log-likelihood as decribed here:
    http://iopscience.iop.org/article/10.1088/0004-637X/750/1/3/pdf
    where the model to fit to data is given by norm*fore_map+const.

    Parameters
    ----------
    norm_guess : float
          initial guess for normalization parameter
    const_guess : float
          initial guess for constant parameter
    fore_map : numpy array
          helapix map of foreground model
    data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See
          next parameter.
    exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in
          counts (beacause foreground map is in flux units by default and it
          needs to be turned to counts to be fitted). While, If data map is
          in flux units, do not declare this parameter, which is None by
          default.
    sr : float or None
          pixel area -> 4*pi/Npix
          
    Returns
    -------
    float
        likelihood value.
        
    """
    a = norm_guess
    b = const_guess
    factorial_data = factorial(data_map)
    lh = 0
    if exp is not None:
        for i, f in enumerate(fore_map):
            lh += (a*f+b)*exp[i]*sr+np.log(factorial_data[i])-data_map[i]*np.log((a*f+b)*exp[i]*sr)
    else:
        for i, f in enumerate(fore_map):
            lh += np.sum(((a*f+b)+np.log(factorial_data[i])-data_map[i]*np.log((a*f+b))))
    return lh


def get_2params_profile_likelihood(lh_matrix, param1_list, param2_list):
    """
    Returns splines with profile likelihood for the two parameters of the fit.
    NOTE: param1 is supposed to be the normalization, param2 the constant.
    """
    
    n_lh = np.amin(lh_matrix, axis=1)
    c_lh = np.amin(lh_matrix, axis=0)
        
    fmt1 = dict(xname='N', xunits='', yname='Likelihood', yunits='')
    spline1 = xInterpolatedUnivariateSplineLinear(param1_list, n_lh, **fmt1)
    fmt2 = dict(xname='C', xunits='', yname='Likelihood', yunits='')
    spline2 = xInterpolatedUnivariateSplineLinear(param2_list, c_lh, **fmt2)
    
    return spline1, spline2
    
def get_param_error(profilelh_spline, param_array, lh_delta=2.3):
    """
    """

    lh_array = profilelh_spline(param_array)
    lh_min_idx = np.argmin(lh_array)
    lherr = lh_array[lh_min_idx]+lh_delta

    if lh_min_idx == 0:
        logger.info('ATT: UPPER limit!')
        sx_err = param_array[0]
        dx_err = param_array[-1]
    elif lh_min_idx == len(lh_array)-1:
        logger.info('ATT: LOWER limit!')
        sx_err = param_array[0]
        dx_err = param_array[-1]
    else:
        sx_err = param_array[np.abs(lh_array[:lh_min_idx]-lherr).argmin()]
        dx_err = param_array[lh_min_idx + np.abs(lh_array[lh_min_idx:]-lherr).argmin()]
        
    return sx_err, dx_err


def fit_foreground_poisson(fore_map, data_map, mask_map=None, n_guess=1.,
                           c_guess=0.1,exp=None, smooth=False, show=False):
    """
    Performs the poissonian fit, recursively computing the log likelihood
    (using poisson_likelihood) for a grid of values of fit parameters around
    the guess. Returns the values of parameters which minimize the log
    likelihood, togather to the 1-sigma error

    Parameters
    ----------
    n_guess : float
          initial guess for normalization parameter
    c_guess : float
          initial guess for constant parameter
    fore_map : numpy array
          helapix map of foreground model
    data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See
          next parameter.
    exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in
          counts (beacause foreground map is in flux units by default and it
          needs to be turned to counts to be fitted). While, If data map is
          in flux units, do not declare this parameter, which is None by
          default.
    smooth : bool
          not implemented yet...
    show : bool
          if true it shows some usefull plot to check if the fit is functioning
          
    Returns
    -------
    float, float, float, float, float, float
        In order: best fit N, best fit C, N's right error, N's left error, 
        C's right error, C's left error
    """
    #show=True
    #mask_map=None
    logger.info('Performing poissonian fit...')
    norm_guess = n_guess
    igrb_guess = c_guess
    nside_out = 64
    mask = 0.
    logger.info('N guess = %.2f - C guess = %.1e'%(norm_guess, igrb_guess))
    if mask_map is None:
        logger.info('fit outside default mask: 30deg gp, 2 deg srcs.')
        mask_f = os.path.join(X_OUT, 'fits/Mask_hp64_src2_gp30.fits')
        mask = hp.read_map(mask_f)
    else:
        logger.info('fit outside mask given in config file.')
        mask = mask_map
        mask = np.array(hp.ud_grade(mask, nside_out=nside_out,
                                      power=-2))
        mask[np.where(mask!=np.amax(mask))[0]] = 0
        mask[np.where(mask==np.amax(mask))[0]] = 1
    logger.info('down grade...')
    fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    data_repix = np.array(hp.ud_grade(data_map, nside_out=nside_out, power=-2))
    _unmask = np.where(mask > 1e-30)[0]

    norm_list = np.linspace(norm_guess-0.8, norm_guess+0.8, 200)
    igrb_list = np.logspace(np.log10(igrb_guess*0.01), np.log10(igrb_guess*10), 200)
    
    logger.info('-------------------------------')
    logger.info('Minimization likelihood run1...')
    lh_list = []
    combinations = list(product(norm_list, igrb_list))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask],
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask], data_repix[_unmask])
            lh_list.append(lh)
    
    lh_list = np.array(lh_list)
    lh_matrix = lh_list.reshape(len(norm_list), len(igrb_list))
    prof_lh_norm, prof_lh_igrb = get_2params_profile_likelihood(lh_matrix, norm_list, igrb_list)
    
    nn = np.linspace(np.amin(norm_list), np.amax(norm_list), 1000)
    cc = np.linspace(np.amin(igrb_list), np.amax(igrb_list), 1000)
    
    lh_min = np.amin(prof_lh_norm(nn))
    logger.info('Minimum -LogL = %s'%lh_min)
    
    norm_min = nn[np.argmin(prof_lh_norm(nn))]
    igrb_min = cc[np.argmin(prof_lh_igrb(cc))]
    logger.info('Run1 results: n=%.3f c=%e'%(norm_min, igrb_min))
    
    norm_sxerr, norm_dxerr = get_param_error(prof_lh_norm, nn, lh_delta=2.3)
    logger.info('Norm err: %.4f - %.4f'%(norm_sxerr, norm_dxerr))
    igrb_sxerr, igrb_dxerr = get_param_error(prof_lh_igrb, cc, lh_delta=2.3)
    logger.info('Igrb err: %.2e - %.2e'%(igrb_sxerr, igrb_dxerr))
    

    """
    logger.info('-------------------------------')
    logger.info('Minimization likelihood run2...')
    norm_list = np.linspace(norm_min-0.3, norm_min+0.3, 100)
    igrb_list = np.linspace(igrb_min*0.1, igrb_min*10, 200)
    lh_list = []
    combinations = np.array(list(product(norm_list, igrb_list)))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask],
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask])
            lh_list.append(lh)
            
    lh_list = np.array(lh_list)
    lh_matrix = lh_list.reshape(len(norm_list), len(igrb_list))
    prof_lh_norm, prof_lh_igrb = get_2params_profile_likelihood(lh_matrix, norm_list, igrb_list)
    
    nn = np.linspace(np.amin(norm_list), np.amax(norm_list), 500)
    cc = np.linspace(np.amin(igrb_list), np.amax(igrb_list), 1000)
    
    lh_min = np.amin(prof_lh_norm(nn))
    lh_delta = lh_min+2.3
    logger.info('Minimum -LogL = %s'%lh_min)
    
    norm_min = nn[np.argmin(prof_lh_norm(nn))]
    igrb_min = cc[np.argmin(prof_lh_igrb(cc))]
    logger.info('Run2 results: n=%.3f c=%e'%(norm_min, igrb_min))
    
    norm_sxerr, norm_dxerr = get_param_error(prof_lh_norm, nn, lh_delta)
    logger.info('Norm err: %.4f - %.4f'%(norm_sxerr, norm_dxerr))
    igrb_sxerr, igrb_dxerr = get_param_error(prof_lh_igrb, cc, lh_delta)
    logger.info('Norm err: %.4f - %.4f'%(igrb_sxerr, igrb_dxerr))
    """
    if show == True:
        
        plt.figure(facecolor='white')
        plt.plot(nn, prof_lh_norm(nn), '-', color='black')
        plt.plot([norm_min, norm_min], [lh_min-10, lh_min+40], color='red')
        plt.plot([norm_sxerr, norm_sxerr], [lh_min-2, lh_min+40], 'r--', alpha=0.7)
        plt.plot([norm_dxerr, norm_dxerr], [lh_min-2, lh_min+40], 'r--', alpha=0.7)
        plt.xlabel('Normalization')
        plt.ylabel('-Log(Likelihood)')
        plt.ylim(lh_min-5, lh_min+30)
        plt.xlim(norm_min-0.2, norm_min+0.2)
        
        plt.figure(facecolor='white')
        plt.plot(cc, prof_lh_igrb(cc), '-', color='black')
        plt.plot([igrb_min, igrb_min], [lh_min-10, lh_min+40], color='red')
        plt.plot([igrb_sxerr, igrb_sxerr], [lh_min-2, lh_min+40], 'r--', alpha=0.7)
        plt.plot([igrb_dxerr, igrb_dxerr], [lh_min-2, lh_min+40], 'r--', alpha=0.7)
        plt.xlabel('Constant')
        plt.ylabel('-Log(Likelihood)')
        plt.ylim(lh_min-5, lh_min+30)
        plt.xlim(igrb_min*0.9, igrb_min*1.1)
        plt.xscale('log')
        
        """
        fig = plt.figure(facecolor='white')
        ax = fig.add_subplot(111)
        
        x, y = np.mgrid(norm_list, igrb_list)
        X, Y = np.mgrid(nn, cc)
        print('---------------', lh_matrix.shape, X.shape, Y.shape)
        print('---------------', lh_matrix.shape, x.shape, y.shape)
        Z = griddata((x, y), lh_matrix, (X, Y), method='linear')


        contours = plt.contour(X, Y, Z, 20, colors='0.4')
        cax = ax.matshow(Z, origin='lower', cmap='RdGy',
                         extent=[np.amin(norm_list), np.amax(norm_list), 
                                np.amin(igrb_list), np.amax(igrb_list)], 
                         aspect='auto', alpha=0.5)
        plt.clabel(contours, inline=True, fontsize=8)
        plt.ylabel('C [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
        plt.xlabel('N')
        ax.xaxis.set_ticks_position('bottom')
        plt.grid('off')
        cb = plt.colorbar(cax, format='$%.1e$')
        cb.set_label('-Log(Likelihood)', rotation=90)
        """
        plt.show()
        
    return norm_min, igrb_min, norm_sxerr, norm_dxerr, igrb_sxerr, igrb_dxerr


def main():
    """Test smodule.
    
    """
    logger.info('No test module is available at the moment... bye bye!')
    return 0









if __name__ == '__main__':
    main()
