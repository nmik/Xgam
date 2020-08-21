#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


"""Analysis module                                                              
"""


import os
import re
import ast
import sys
import argparse
import numpy as np
import healpy as hp
from itertools import combinations
from scipy.optimize import curve_fit

__description__ = 'Makes the analysis'


"""Command-line switches.                                                       
"""


from Xgam import X_OUT, X_CONFIG
from Xgam.utils.PolSpice_ import *
from Xgam.utils.wbeamfunc_ import *
from Xgam.utils.logging_ import logger

EMIN_STR = '(?<=\_)[0-9]+(-)'
EMAX_STR = '(?<=\-)[0-9]+(\.fits)'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--show', type=ast.literal_eval, choices=[True, False],
                    default=False, help='True if you want to see the maps')
PARSER.add_argument('--overwrite', type=bool, choices=[True, False],
                    default=False, help='if True overwrite existing output files')

if (sys.version_info > (3, 0)):
	from importlib.machinery import SourceFileLoader
	def get_var_from_file(filename):
		f = open(filename)
		global data
		data = SourceFileLoader('data', filename).load_module()
		f.close()
else:
	import imp
	def get_var_from_file(filename):
		f = open(filename)
		global data
		data = imp.load_source('data', '', f)
		f.close()

def mkAuto(**kwargs):
    """                                      
    """
    get_var_from_file(kwargs['config'])

    maps_ = data.MAPS_LIST
    masks_ = data.MASKS_LIST
    fermi_wb_matrix = data.FERMI_WBEAM_MATRIX
    lss_wb_ = data.LSS_TRACER_WBEAM_LIST
    out_label = data.OUT_LABEL
    gamma = data.GAMMA
    l_max_ = data.MAX_APS_MULTIPOLE_FIT
    lmax = data.BINNING_MAX_MULTIPOLE
    lmin = data.BINNING_MIN_MULTIPOLE
    bin_num = data.BINNING_MULTIPOLE_NBIN
    bin_alg = data.BINNING_MULTIPOLE_ALGORITHM
    bin_custom = data.BINNING_CUSTOM
    fit_function = data.FIT_FUNCTION
    cn_list = data.FERMI_CN_LIST
    
    if type(maps_) == str and maps_.endswith(('.txt','.dat')):
        maps_ = open(maps_, 'r')
        maps_ = maps_.read().splitlines()
    
    if type(masks_) == str and masks_.endswith(('.txt','.dat')):
        masks_ = open(masks_, 'r')
        masks_ = masks_.read().splitlines()

    logger.info('Starting ...')
    cl_txt_f = os.path.join(X_OUT,'%s_autocorrelation.txt'%(out_label))
    if os.path.exists(cl_txt_f):
        logger.info('ATT: Output file already exists!')
        logger.info(cl_txt_f)
        if not kwargs['overwrite']:
            sys.exit()
        else:
            logger.info('ATT: Overwriting files!')
            
    cl_txt = open(cl_txt_f, 'w')
    wl_ = []
    logger.info('----- AUTO -----')
    for i, m1_f in enumerate(maps_):
        l_max = l_max_[i]
        logger.info('Considering map: %s'%m1_f)
        m1 = hp.read_map(m1_f)
        
        if len(masks_) == len(maps_):
            mask_f = masks_[i]
        elif len(masks_) == 1:
            mask_f = masks_[0]
        else:
            logger.info('ERROR: Inconsistent number of masks!'+
                         'It should be either 1 or the same number as Fermi maps!')
            sys.exit()
        logger.info('Considering Fermi mask: %s'%mask_f)
        mask = hp.read_map(mask_f)
        
        m1_masked = hp.ma(m1)
        m1_masked.mask = np.logical_not(mask)
        
        logger.info('Computing the pixel window function ...')
 
        nside = hp.npix2nside(len(m1))
        wpix = hp.wpix = hp.sphtfunc.pixwin(nside, lmax=lmax-1)
        cn = cn_list[i]
        fsky = np.sum(mask)/len(mask)
		
        if fermi_wb_matrix is not None:
            logger.info('Computing the fermi beam window function ...')
            spec = get_powerlaw_spline(gamma)
            emin_str = re.search(r'%s'%EMIN_STR, m1_f)
            emin = float(emin_str.group(0).replace('-',''))
            emax_str = re.search(r'%s'%EMAX_STR, m1_f)
            emax = float(emax_str.group(0).replace('.fits', ''))
            wb_spline = get_1D_wbeam(fermi_wb_matrix, spec, e_min=emin, e_max=emax)
            wb = wb_spline(np.arange(lmax))
            
            if kwargs['show']:
                wb_spline.plot(show=False, label='%.2f-%.2f GeV'%(emin/1000, emax/1000), color='0.3')
                plt.ylabel('W$_{beam}$', size=15)
                plt.xlabel('$l$', size=15)
                plt.legend(loc=1, fontsize=14)
                plt.grid('off')
            
        # ------- TO BE DEFINED ------
        if lss_wb_ is not None:
            logger.info('Getting the LSS beam window function ...')
            lss_wb = np.ones(len(wb))
            
        cl_txt.write('# AUTOCORR ---> %i x %i\n'%(i, i))

        if lss_wb_ is None:
            wl = wb * wpix
            wl_.append(wl)
        else:
            # ------- TO BE DEFINED ------
            wl = lss_wb * wpix
            wl_.append(wl)
            
        out_folder =  os.path.join(X_OUT, 'output_pol')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        pol_dict = data.POLCEPICE_DICT
        for key in pol_dict:
            if key == 'clfile':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_cl.txt'%(out_label, i))
            if key == 'cl_outmap_file':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_clraw.txt'%(out_label, i))
#             if key == 'covfileout':
#                 pol_dict[key] = os.path.join(out_folder,'%s_%i_cov.fits'%(out_label, i))
            if key == 'mapfile':
                pol_dict[key] = m1_f
            if key == 'maskfile':
                pol_dict[key] = mask_f
            if key == 'mapfile2':
                pol_dict[key] = 'NO'
            if key == 'maskfile2':
                pol_dict[key] = 'NO'
        config_file_name = os.path.join(out_folder, 'pol_autoconfig_%s_%i.txt'%(out_label, i))
       
        _l, _cl, _cl_err = pol_cl_calculation(pol_dict, config_file_name, cov=False,
                                                        wl_array=wl, cn=cn, rebin=True,
                                                        nbin=bin_num, bin_type=bin_alg,
                                                        lmin=lmin, lmax=lmax,
                                                        custom_bins=bin_custom)
                                                                
        fit_idx = np.where((_l >= 30)&(_l <= l_max))[0]

        try: 
            popt, pcov = curve_fit(fit_function, _l[fit_idx], _cl[fit_idx], sigma=_cl_err[fit_idx])
            logger.info('\n')
            logger.info('********************************')
            logger.info('Multipole fit range: %i - %i'%(30, l_max))
            logger.info('Best fit param: %s'%str(popt))
            logger.info('********************************\n')
            cl_txt.write('Cn : %s\n'%str(cn))
            cl_txt.write('Multipole fit range: %i - %i\n'%(30, l_max))
            cl_txt.write('Best fit params : %s\n'%str(popt))
             
            
            plt.figure()
            plt.errorbar(_l[fit_idx], _cl[fit_idx], yerr=_cl_err[fit_idx], fmt='.', 
                                            label='%.2f-%.2f GeV'%(emin/1000, emax/1000))
            plt.plot(_l[fit_idx], fit_function(_l[fit_idx], *popt))
            if bin_alg == 'log':
                plt.xscale('log')
            plt.xlabel('$\ell$')
            plt.ylabel('$C_{\ell}$')
            plt.legend()
            plt.savefig('output/figs/Autocorr_%ix%i'%(i,i))
                 
            if kwargs['show']:  
                plt.show()
             
        except RuntimeError:
            logger.info('\n')
            logger.info('********************************')
            logger.info("Error - curve_fit failed")
            logger.info('********************************\n')
            cl_txt.write('Cn : %s\n'%str(cn))
            cl_txt.write('Error - curve_fit failed\n')

#             plt.figure()
#             plt.errorbar(_l[fit_idx], _cl[fit_idx], yerr=_cl_err[fit_idx], fmt='.')
#             if bin_alg == 'log':
#                 plt.xscale('log')    
#             plt.show()
    
    
    
    
    logger.info('----- CROSS -----')
    
    for i, m1_f in enumerate(maps_):
        wl1 = wl_[i]
        mask1_f = masks_[i]
        for j in range(i, len(maps_)):
            wl2 = wl_[j]
            m2_f = maps_[j]
            mask2_f = masks_[j]
            if m1_f == m2_f:
               continue
            
            out_folder =  os.path.join(X_OUT, 'output_pol')
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            pol_dict = data.POLCEPICE_DICT
            for key in pol_dict:
                if key == 'clfile':
                    pol_dict[key] = os.path.join(out_folder,'%s_%i%i_cl.txt'%(out_label, i, j))
                if key == 'cl_outmap_file':
                    pol_dict[key] = os.path.join(out_folder,'%s_%i%i_clraw.txt'%(out_label, i, j))
                if key == 'mapfile':
                    pol_dict[key] = m1_f
                if key == 'maskfile':
                    pol_dict[key] = mask1_f
                if key == 'mapfile2':
                    pol_dict[key] = m2_f
                if key == 'maskfile2':
                    pol_dict[key] = mask2_f
            config_file_name = os.path.join(out_folder, 'pol_autocrossconfig_%s_%i%i.txt'%(out_label, i, j))
       
            wl = np.sqrt(wl1*wl2)
            _l, _cl, _cl_err = pol_cl_calculation(pol_dict, config_file_name, cov=False,
                                                        wl_array=wl, cn=0, rebin=True,
                                                        nbin=bin_num, bin_type=bin_alg,
                                                        lmin=lmin, lmax=lmax,
                                                        custom_bins=bin_custom)
            
            cl_txt.write('# CROSSCORR ---> %i x %i\n'%(i, j))
            
            fit_idx = np.where((_l >= 30)&(_l <= l_max_[i]))[0]

            try: 
                popt, pcov = curve_fit(fit_function, _l[fit_idx], _cl[fit_idx], sigma=_cl_err[fit_idx])
                logger.info('\n')
                logger.info('********************************')
                logger.info('Multipole fit range: %i - %i'%(30, l_max_[i]))
                logger.info('Best fit param: %s'%str(popt))
                logger.info('********************************\n')
                cl_txt.write('Cn : 0')
                cl_txt.write('Multipole fit range: %i - %i\n'%(30, l_max_[i]))
                cl_txt.write('Best fit params : %s\n'%str(popt))
             
            
                plt.figure()
                plt.errorbar(_l[fit_idx], _cl[fit_idx], yerr=_cl_err[fit_idx], fmt='.', 
                                            label='%.2f-%.2f GeV'%(emin/1000, emax/1000))
                plt.plot(_l[fit_idx], fit_function(_l[fit_idx], *popt))
                if bin_alg == 'log':
                    plt.xscale('log')
                plt.xlabel('$\ell$')
                plt.ylabel('$C_{\ell}$')
                plt.legend()
                plt.savefig('output/figs/Autocorr_%ix%i'%(i,j))
                 
                if kwargs['show']:  
                    plt.show()
             
            except RuntimeError:
                logger.info('\n')
                logger.info('********************************')
                logger.info("Error - curve_fit failed")
                logger.info('********************************\n')
                cl_txt.write('Error - curve_fit failed\n')

            
    cl_txt.close()

    logger.info('Created %s'%cl_txt_f)



if __name__ == '__main__':
    args = PARSER.parse_args()
    mkAuto(**args.__dict__)
