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


"""general parsing functions.
"""

import os
import numpy as np
from astropy.io import fits as pf
from matplotlib import pyplot as plt

from Xgam import X_CONFIG
from Xgam.utils.logging_ import logger
from Xgam.utils.matplotlib_ import *
from Xgam.utils.spline_ import xInterpolatedBivariateSplineLinear
from Xgam.utils.spline_ import xInterpolatedUnivariateSplineLinear


def parse_polspice_aps(cl_file):
    """
    Parsing of the *_cls.txt files.

    cl_file : str
        This file is created by bin/mkpolspicecl.py app.
    """
    ls = []
    cls = []
    clerrs = []
    _cov_ = []
    emin, emax = [], []
    f = open(cl_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2 = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
        if 'multipole\t' in line:
            l = np.array([float(item) for item in line.split()[1:]])
            ls.append(l)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
        if 'Cl_ERR' in line:
            cl_err = np.array([float(item) for item in line.split()[1:]])
            clerrs.append(cl_err)
        if 'COV_FILE' in line:
            cov_f = line.split()[-1]
            try:
                cov = np.load(cov_f, allow_pickle=True)
                _cov_.append(cov)
            except:
                pass
    f.close()
    return [np.array(emin),np.array(emax), np.array(ls), np.array(cls), \
                                                  np.array(clerrs), np.array(_cov_)]

def parse_datafluxmaps(out_dataflux_map_file):
    """Parsing of *_parameters.txt files.

       cl_param_file : str
           This file is created by bin/mkdatarestyle.py app.
    """
    logger.info('Loading parameters from %s'%out_dataflux_map_file)
    ff = open(out_dataflux_map_file, 'r')
    _emin, _emax, _emean = [], [], []
    _f, _ferr, _cn, _fsky = [], [], [], []
    _n, _n_sx, _n_dx, _c, _c_sx, _c_dx = [], [], [], [], [], []
    for line in ff:
        try:
            emin, emax, emean, f, ferr, cn, fsky, n, nsx, ndx, c, csx, cdx = [float(item) for item in \
                                                        line.split()]  
            _emin.append(emin)
            _emax.append(emax)
            _emean.append(emean)
            _f.append(f)
            _ferr.append(ferr)
            _cn.append(cn)
            _fsky.append(fsky)
            _n.append(n)
            _n_sx.append(nsx)
            _n_dx.append(ndx)
            _c.append(c)
            _c_sx.append(csx)
            _c_dx.append(cdx)
        except:
            try:
                emin, emax, emean, f, ferr, cn, fsky = [float(item) for item in line.split()]   
                _emin.append(emin)
                _emax.append(emax)
                _emean.append(emean)
                _f.append(f)
                _ferr.append(ferr)
                _cn.append(cn)
                _fsky.append(fsky)   
            except:
                pass
            
    ff.close()
    return [np.array(_emin), np.array(_emax), np.array(_emean), np.array(_f), \
        np.array(_ferr), np.array(_cn), np.array(_fsky), np.array(_n),np.array(_n_sx),\
        np.array(_n_dx), np.array(_c), np.array(_c_sx), np.array(_c_dx)]

def get_energy_from_fits(fits_file, minbinnum=0, maxbinnum=100, mean='log'):
    """Returns a list with the center values of the energy bins

       fits_file: str
           fits file, usually we want to do to this at the level of
           gtbin output file
       mean: str
           'log' or 'lin', depending on the algorithm to use to 
           compute the mean
    """
    f = pf.open(fits_file)
    ebounds = f[1].data
    _emin = ebounds['E_MIN'][minbinnum:maxbinnum+1]/1000
    _emax = ebounds['E_MAX'][minbinnum:maxbinnum+1]/1000
    emean = []        
    if mean == 'log':
        for emin, emax in zip(_emin, _emax):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(_emin, _emax):
            emean.append(0.5*(emin+emax))
    f.close()
    return np.array(_emin), np.array(_emax), np.array(emean)
    
    
def get_psf_th_en_bivariatespline(PSF_FILE, show=False):
    """ Get the PSF from the fits file created by gtpsf       
    
        PSF_FILE: str                                                                                                                                                     
           file .fits generated by gtpsf                                                                                                           
    """
    hdu_list = pf.open(PSF_FILE)
    _th = np.radians(hdu_list['THETA'].data.field('Theta'))
    _en = hdu_list['PSF'].data.field('ENERGY')
    _psf = hdu_list['PSF'].data.field('PSF')
    fmt = dict(yname='Theta', yunits='rad', xname='Energy',
               xunits='MeV', zname='PSF')
    psf_th_e_spline = xInterpolatedBivariateSplineLinear(_en, _th, _psf, **fmt)
    hdu_list.close()
    if show == True:
    	plt.imshow(np.log10(_psf))
    	plt.figure()
    	for e in _en[::50]:
    		psf_th = psf_th_e_spline.vslice(e)
    		fmt = dict(xname='cos(th)', xunits='', yname='psf[cos(th)]', yunits='')
    		plt.plot(np.degrees(psf_th.x), psf_th.y, label='%.2f'%e)
    	plt.yscale('log')
    	plt.legend()
    	plt.figure()
    	psf_e = np.sum(_psf, axis=1)
    	fmt = dict(xname='cos(th)', xunits='', yname='psf[cos(th)]', yunits='')
    	plt.plot(_en, 1/psf_e)
    	plt.yscale('log')
    	plt.xscale('log')
    	plt.show()
    return psf_th_e_spline
    
    
def get_psf_en_univariatespline(PSF_FILE, show=False):
	""" Gives the 95% containment angle as a function of energy.
	
	    PSF_FILE: str
	       PSF file produced with gtpsf.
	"""
	hdu_list = pf.open(PSF_FILE)
	_th = np.radians(hdu_list['THETA'].data.field('Theta'))
	_en = hdu_list['PSF'].data.field('ENERGY')
	_psf = hdu_list['PSF'].data.field('PSF')
	_psf_e = []
	for j, en in enumerate(_en):
		psf_ph_spline = xInterpolatedUnivariateSplineLinear(_th, _psf[j])
		tot_integral = psf_ph_spline.integral(_th[0], _th[-1])
		k = 1
		while psf_ph_spline.integral(_th[0], _th[k]) < 0.95*tot_integral:
			k = k+1
		_psf_e.append(np.degrees(_th[k]))
	_psf_e = np.array(_psf_e)
	psf_e_spline = xInterpolatedUnivariateSplineLinear(_en,  _psf_e)
	if show == True:
		psf_e_spline.plot(show=False)
		plt.xscale('log')
		plt.show()
	return psf_e_spline

	

if __name__ == '__main__': 
	""" test module
	"""
	ebin_f = 'output/ebinning.fits'
	emin, emax, emean = get_energy_from_fits(ebin_f)
	print(emean)
	
	psf_f = os.path.join(X_CONFIG, 'fits/psf_SV_t32.fits')
	psf_bispline = get_psf_th_en_bivariatespline(psf_f, show=True)
	psf_unispline = get_psf_en_univariatespline(psf_f, show=True)
