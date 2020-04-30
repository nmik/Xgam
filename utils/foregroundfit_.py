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

import os
import re
import numpy as np
import healpy as hp
from numba import jit
import astropy.io.fits as pf
from itertools import product
from scipy.special import factorial
import matplotlib.pyplot as plt

from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg

FORE_EN = re.compile('\_\d+\.')

def get_fore_integral_flux_map(fore_files_list, e_min, e_max):
    """Returns the foreground map integrated between e_min and e_max
       A powerlaw is assumed fore the foreground energy spectrum, hence
       the interpolation between 2 given maps at given energies (given
       by the model) is done in logarithmic scales.

       fore_files_list: list of str
           Ordered list of the foreground files (one for each energy)
       e_min: float
           the min of the energy bin
       e_max: float
           the max of the energy bin
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
    """Returns the first element on the right and the first on the left
       of a given value (en_val), among all values in an ordered array
       (en_arr).

       en_val : float
           mean energy
       en_arr : float
           array of the energies at which the foreground model is given.
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
def poisson_likelihood(norm_guess, const_guess, fore_map, data_map, exp=None, sr=None):
    """Compute the log-likelihood as decribed here:
       http://iopscience.iop.org/article/10.1088/0004-637X/750/1/3/pdf
       where the model to fit to data is given by norm*fore_map+const.

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
    """
    a = norm_guess
    b = const_guess
    factorial_data = factorial(data_map)
    lh = 0
    if exp is not None:
        for i, f in enumerate(fore_map):
            lh += (a*f+b)*exp[i]*sr + np.log(factorial_data[i]) - \
                data_map[i]*np.log((a*f+b)*exp[i]*sr)
    else:
        for i, f in enumerate(fore_map):
            lh += np.sum(((a*f+b)+np.log(factorial_data[i]) -\
                              data_map[i]*np.log((a*f+b))))
    return lh


def fit_foreground_poisson(fore_map, data_map, mask_map=None, n_guess=1.,
                           c_guess=0.1,exp=None, smooth=False, show=False):
    """Performs the poisonian fit, recursively computing the log likelihood
       (using poisson_likelihood) for a grid of values of fit parameters around
       the guess. Returns the values of parameters which minimize the log
       likelihood, togather to the 1-sigma error

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
    """
    #show=True
    logger.info('Performing poissonian fit...')
    norm_guess = n_guess
    igrb_guess = c_guess
    nside_out = 64
    mask = 0.
    if mask_map is None:
        logger.info('fit outside default mask: 30deg gp, 2 deg srcs.')
        mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask64_src2_gp30.fits')
        mask = hp.read_map(mask_f)
    else:
        logger.info('fit outside mask given in config file.')
        mask = mask_map
    logger.info('down grade...')
    fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    data_repix = np.array(hp.ud_grade(data_map, nside_out=nside_out,
                                      power=-2))
    mask_repix = np.array(hp.ud_grade(mask, nside_out=nside_out,
                                      power=-2))
    mask_repix[np.where(mask_repix!=np.amax(mask_repix))[0]] = 0
    mask_repix[np.where(mask_repix==np.amax(mask_repix))[0]] = 1
    _unmask = np.where(mask_repix > 1e-30)[0]
    norm_list = np.linspace(norm_guess*0.3, norm_guess*1.5, 50)
    igrb_list = np.linspace(igrb_guess*0.01, igrb_guess*20., 200)
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
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask])
            lh_list.append(lh)
    lh_min = np.argmin(np.array(lh_list))
    (norm_min, igrb_min) = combinations[lh_min]
    logger.info('Run1 results: n=%.3f c=%.1e'%(norm_min, igrb_min))
    norm_list = np.linspace(norm_min*0.7, norm_min*1.2, 51)
    igrb_list = np.linspace(igrb_min*0.1, igrb_min*10, 101)
    logger.info('Minimization likelihood run2...')
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
    lh_min = np.argmin(lh_list)
    logger.info('Minimum -LogL = %s'%lh_list[lh_min])
    (norm_min, igrb_min) = combinations[lh_min]
    logger.info('Run2 results: n=%.3f c=%e'%(norm_min, igrb_min))
    lh_delta = np.array(lh_list)[lh_min]+2.3
    index = np.where(np.array(lh_list) < lh_delta)[0]
    _norm = np.array([x[0] for x in combinations[index]])
    logger.info('Norm err: %.4f - %.4f'%(_norm[0], _norm[-1]))
    _igrb = np.array([x[1] for x in combinations[index]])
    logger.info('Igrb err: %.3e - %.3e'%(np.amin(_igrb), np.amax(_igrb)))

    if show == True:
        n = np.array([x[0] for x in combinations])
        plt.figure(facecolor='white')
        plt.plot(n, lh_list, 'o', color='coral', alpha=0.3)
        plt.plot(norm_min, lh_list[lh_min] , 'r*')
        plt.plot([_norm[0], _norm[-1]], [lh_delta, lh_delta], 'r-')
        plt.xlabel('Normalization')
        plt.ylabel('-Log(Likelihood)')

        igrb = np.array([x[1] for x in combinations])
        plt.figure(facecolor='white')
        plt.plot(igrb, lh_list, 'o', color='coral', alpha=0.3)
        plt.plot(igrb_min, lh_list[lh_min] , 'r*')
        plt.plot([np.amin(_igrb), np.amax(_igrb)], [lh_delta, lh_delta], 'r-')
        plt.xlabel('Constant')
        plt.ylabel('-Log(Likelihood)')

        fig = plt.figure(facecolor='white')
        z = lh_list
        zmin = lh_list[lh_min]
        z.shape = (len(norm_list), len(igrb_list))
        ax = fig.add_subplot(111)
        cax = ax.matshow(z, origin='lower', cmap='Spectral',
                    aspect='auto')
        plt.xlabel('C [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
        plt.ylabel('N')
        x_ticks = np.linspace(np.amin(igrb_list), np.amax(igrb_list), 6)
        formatting_function = np.vectorize(lambda f: format(f, '6.1E'))
        x_ticks = list(formatting_function(x_ticks))
        y_ticks = list(np.around(np.linspace(np.amin(norm_list),
                                             np.amax(norm_list), 6),
                                 decimals=3))
        ax.set_yticklabels(['']+y_ticks)
        ax.set_xticklabels(['']+x_ticks)
        ax.xaxis.set_ticks_position('bottom')
        plt.grid('off')
        cb = plt.colorbar(cax, format='$%.1e$')
        cb.set_label('-Log(Likelihood)', rotation=90)
        norm_min_ind = list(norm_list).index(norm_min)
        igrb_min_ind = list(igrb_list).index(igrb_min)
        _norm_ind = []
        _igrb_ind = []
        for i in range(0, len(index)):
            _norm_ind.append(list(norm_list).index(_norm[i]))
            _igrb_ind.append(list(igrb_list).index(_igrb[i]))
        _norm_ind = np.array(_norm_ind)
        _igrb_ind = np.array(_igrb_ind)
        plt.contourf(z, [zmin, zmin+2.3, zmin+4.61, zmin+5.99],
                     colors='w', origin='lower', alpha=0.3)
        plt.show()
    return norm_min, igrb_min, _norm[0], _norm[-1], np.amin(_igrb), \
        np.amax(_igrb)


def main():
    """Test session
    """
    logger.info('No test module is available at the moment... bye bye!')
    return 0









if __name__ == '__main__':
    main()
