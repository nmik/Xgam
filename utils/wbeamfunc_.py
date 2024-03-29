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


"""Functions to generate Wbeam(l, E) matrix.
"""

import os
import numpy as np
import scipy.special as sp
import astropy.io.fits as pf

from Xgam import X_OUT
from Xgam.utils.logging_ import logger
import matplotlib.pyplot as plt 
from Xgam.utils.spline_ import xInterpolatedBivariateSplineLinear
from Xgam.utils.spline_ import xInterpolatedUnivariateSplineLinear
from Xgam.utils.parsing_ import get_energy_from_fits

def get_pl_vs_th(l_, costh_):
    """ 
    Gives Pl values correspoding to a given _costh array
    
    Parameters
    ----------
    l_ : array
        array of multipoles ell (integers from 0 to ell_max)
    costh_ : array
        array of cosine theta
       
    Returns
    -------
    array
         Pl values correspoding to a given _costh array
         
    """
    pl_th_ = sp.eval_legendre(l_, costh_)
    
    return pl_th_
    

def build_wpix(nside, l_max=1500):
    """
    Calculate the pixel window function.
    
    Parameters
    ----------
    nside : int
 	    nside of the map you are analyzing
    l_max : int
        the maximum multipol at wich to calculate the pixel window function
           
    Returns
    -------
    spline
        spline of the pixel window function  
                
    """
    wpix = hp.sphtfunc.pixwin(nside)[:l_max]
    fmt = dict(xname='$l$', xunits='', yname='W_{pixel}', yunits='')
    wpix_spline = xInterpolatedUnivariateSplineLinear(np.arange(l_max), wpix, **fmt)
                                                      
    return wpix_spline
    
    
def build_wbeam(psf, l_, out_file):
    """
    Calculates the Wbeam(l, E) and return a bivariate slpine.
    
    Parameters
    ----------
    psf : spline
        psf bivariate spline generated by get_psf
    l_  :  array
        array of multipoles at wicht to compute the Wbeam
    out_file : str
        name of the output txt file that will be created in output/
             
    Returns
    ------- 
    spline
        spline of the beam window function   
        
    """
    out_txt = open(os.path.join(X_OUT, out_file), 'w')
    en_ = psf.x
    _wb_ = []
    energy = str(list(en_)).replace('[','').replace(']','').replace(', ', ' ')
    out_txt.write('l\t%s\n'%energy)
    for e in en_:
        wb_e = np.array([])
        psf_th = psf.vslice(e)
        for l in l_:
            pl_th = get_pl_vs_th(l, np.cos(psf_th.x))
            fmt = dict(xname='th', xunits='rad', yname='convolution', \
                           yunits='MeV')
            _conv = xInterpolatedUnivariateSplineLinear(psf_th.x, \
                                    np.sin(psf_th.x)*psf_th.y*pl_th, **fmt)
            wb_e_l = min(1., 2*np.pi*(_conv.integral(np.amin(psf_th.x), \
                                                         np.amax(psf_th.x))))
            logger.info('Wbeam(%i, %.2f) = %e'%(l, e, wb_e_l))
            wb_e = np.append(wb_e, [wb_e_l])
        _wb_.append(wb_e)
    _wb_ = np.array(_wb_)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                             yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear(l_, en_, _wb_.T, **fmt)
    for i, l in enumerate(l_):
        wb = str(list(_wb_.T[i])).replace('[','').replace(']','').\
            replace(', ', ' ')
        out_txt.write('%i\t%s\n'%(l, wb))
    out_txt.close()
    return _wb_
    

def wbeam_parse(wb_file, l_max=1500):
    """
    Created to parse the txt file given in output by build_wbeam.
    
    Parameters
    ----------
    wb_file : str 
        .txt file generated by build_wbeam
    l_max : int
        the maximum multipol at wich to calculate the beam window function
        
    Returns
    -------
    array, array, tensor
        In order: array of the energies, array of the multipoles, tensor of the 2D Wbeam(l, E)

    """
    f = open(wb_file, 'r')
    e_ = []
    l_ = []
    _z_ = []
    for i, line in enumerate(f):
        if 'l' in line:
            e_.extend(float(item) for item in line.split()[1:])
        elif str(i-1) == line.split()[0]:
            l_.append(i+1)
            _z_.append(np.array([float(item) for item in line.split()[1:]]))
    return np.array(e_), l_, np.array(_z_)
    
    
def get_2D_wbeam(wb_file, show=False):
    """ Retrive the bivariate spline of the Wbeam function recorded in a txt file.
    
    Parameters
    ----------
    wb_file : str 
    	.txt file generated by build_wbeam
    show : bool
    	if True the Wbeam matrix is plotted
             
    Returns
    -------
    spline
        2D spline of the beam window function  
    	
    """
    en_, l_, _z_ = wbeam_parse(wb_file)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                             yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear(l_, en_, _z_, **fmt)
    
    if show == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(_z_.T, origin='lower', aspect='auto', cmap='viridis')
        #en_tick = ['1$\cdot$10$^{2}$', '4.6$\cdot$10$^{2}$', 
        #           '1$\cdot$10$^{3}$','1.$\cdot$10$^{4}$', 
        #           '4.6$\cdot$10$^{4}$', '2.5$\cdot$10$^{5}$']
        #l_tick = ['0','200','400','600','800','1000','1200','1400']
        #ax.set_yticklabels(['']+en_tick, size=13)
        #ax.set_xticklabels(['']+l_tick, size=13)
        ax.xaxis.set_ticks_position('bottom')
        #plt.title('PSF 1+2+3', size=18)
        plt.xlabel('Multipole $l$', size=15)
        plt.ylabel('Energy', size=15)
        cbar = plt.colorbar(cax)
        cbar.set_label('$W_{beam}(l, E)$', size=15)
        plt.grid(False)
        plt.show()
    
    return wbeam
    

def get_powerlaw_spline(index=2.3):
    """
    Returns a simple power law spline with a given spectral index.
    
    Parameters
    ----------
    index : float
    	spectral index of the powerlaw you want to generate
             
    Returns
    -------
    spline
        spline of the energy spectrum (power law with index = -gamma)
    	
    """
    EMIN = 2
    EMAX = 6
    en_ = np.logspace(EMIN, EMAX, 2000)
    pl_ = en_**(-index)
    fmt = dict(xname='$E$', xunits='MeV', yname='E$^{-%.1f}$'%index,
                   yunits='')
    spectrum = xInterpolatedUnivariateSplineLinear(en_, pl_, **fmt)
    
    return spectrum
    
    
def get_1D_wbeam(wb_file, spectrum_spline, e_min, e_max):
    """
    Returns the integral Wbeam as a function of the multiple in a certain energy interval.
       
    Parameters
    ----------
    wb_file : str 
        .txt file generated by build_wbeam
    spectrum_spline : numpy spline 
        spline of the spectrum returned by get_powerlaw_spline
    e_min : float
	    lower limit of the energy interval 
    e_max : float
        upper limit of the energy interval 
       
    Returns
    -------
    spline
        spline Wbeam(l) resulting from Wbeam(E, l) integrated between Emin and Emax
        
    """
    en_, l_, _z_ = wbeam_parse(wb_file)
    int_spec = spectrum_spline.integral(e_min, e_max)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    int_wbeam = xInterpolatedBivariateSplineLinear( l_, en_, _z_*spectrum_spline(en_),
                                                   **fmt)
    int_wbeam_l = []
    for l in l_:
        l_slice = int_wbeam.vslice(l)
        int_wbeam_l.append(l_slice.integral(e_min, e_max))
    int_wbeam_l = np.array(int_wbeam_l)
    
    fmt2 = dict(xname='$l$', xunits='', 
                yname='integral W$_{beam}$ ($%.1f - %.1f$)'%(e_min, e_max), yunits='')
    wbeam = xInterpolatedUnivariateSplineLinear(l_, int_wbeam_l/int_spec, **fmt2)
    norm = 1/wbeam.y[0]
    return wbeam
 
    
def main():
    """Test module
    
    """
    
    CREATE = False
    RETRIEVE = True
     
    # To build a new set of Wbeam functions:
    if CREATE:
        from Xgam.utils.parsing_ import get_psf_th_en_bivariatespline
        
        psf_f = os.path.join(X_OUT, 'fits/psf_SV_t32.fits')
        psf = get_psf_th_en_bivariatespline(psf_f)
        out_wbeam_txt = 'P8R3_SOURCEVETO_V2_evt32_wbeam.txt'
        l_ = np.arange(2100)
        wb = build_wbeam(psf, l_, out_wbeam_txt)
    
    if RETRIEVE:
        import healpy as hp
    
        out_wbeam_txt = os.path.join(X_OUT, 'P8R3_SOURCEVETO_V2_evt56_wbeam.txt')
        wb_2d = get_2D_wbeam(out_wbeam_txt, show=True)
        
        wpix = hp.sphtfunc.pixwin(256)
        
        logger.info('Wait...computing integral Wbeam at some energy intervals!')
        gamma = 2.3
        spec = get_powerlaw_spline(gamma)
        
#         _emin = np.array([631.0,1202.3,2290.9,4786.3,9120.1,17378.0,36307.8,69183.1,131825.7])
#         _emax =  np.array([1202.3,2290.9,4786.3,9120.1,17378.0,36307.8,69183.1,131825.7,1000000])
        _emin = np.array([1000.0, 2089.30, 4786.30, 10964.78])
        _emax =  np.array([2089.30, 4786.30, 10964.78, 25118.86])
        
        c = ['0.1', '0.2','0.4','0.5','0.6','0.7','0.8','0.9','0.95']
        c = ['0.1', '0.4','0.6','0.8','0.95']

        fig = plt.figure(facecolor='white')
        
        wb_ = []
        plot = fig.add_subplot(111)
        plt.plot(np.arange(768), wpix[:768], '--', color='red', label='W$_{pix}$ (NSIDE=256)')
        for i, (emin, emax) in enumerate(zip(_emin, _emax)):
            wb = get_1D_wbeam(out_wbeam_txt, spec, e_min=emin, e_max=emax)
            wb_.append(wb.y)
            print('wl(l=1)**2 = %e'%(wb(2)*wpix[1])**2)
            wb.plot(show=False, label='%.2f-%.2f GeV'%(emin/1000, emax/1000), color='%s'%c[i])
            
        plt.ylabel('W$_{beam}$', size=18)
        plt.xlabel('$l$', size=18)
        plt.ylim(-1,1.1)
        plt.xlim(0, 1500)
        plt.legend(loc=3, fontsize=15, fancybox=True)
        plot.tick_params(axis='both', which='major', labelsize=14)
        #plt.grid('off')
        
        txt_f = 'output/Wbeam_SV_t56_gxn_allbins.txt'
        np.savetxt(txt_f, np.array(wb_).T, header='ell, wbeam(l) (E1, E2, E3, E4)')
#         np.savetxt(txt_f, np.array(wb_).T, header='E1\tE2\tE3\tE4\tE5\tE6\tE7\tE8\tE9')

        
        
        plt.show()
    
if __name__ == '__main__':
    main()
