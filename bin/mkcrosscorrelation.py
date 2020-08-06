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

def mkCross(**kwargs):
    """                                      
    """
    get_var_from_file(kwargs['config'])

    fermi_maps_ = data.FERMI_MAPS_LIST
    fermi_masks_ = data.FERMI_MASKS_LIST
    fermi_wb_matrix = data.FERMI_WBEAM_MATRIX
    lss_maps_ = data.LSS_TRACER_MAPS_LIST
    lss_masks_ = data.LSS_TRACER_MASK_LIST
    lss_wb_ = data.LSS_TRACER_WBEAM_LIST
    out_label = data.OUT_LABEL
    bin_label = data.BINNING_LABEL
    gamma = data.GAMMA
    l_max = data.MAX_APS_MULTIPOLE
    lmax = data.BINNING_MAX_MULTIPOLE
    lmin = data.BINNING_MIN_MULTIPOLE
    bin_num = data.BINNING_MULTIPOLE_NBIN
    bin_alg = data.BINNING_MULTIPOLE_ALGORITHM
    bin_custom = data.BINNING_CUSTOM
    
    if type(fermi_maps_) == str and fermi_maps_.endswith(('.txt','.dat')):
        fermi_maps_ = open(fermi_maps_,'r')
        fermi_maps_ = fermi_maps_.read().splitlines()
    
    if type(fermi_masks_) == str and fermi_masks_.endswith(('.txt','.dat')):
        fermi_masks_ = open(fermi_masks_,'r')
        fermi_masks_ = fermi_masks_.read().splitlines()
    
    if type(lss_maps_) == str and lss_maps_.endswith(('.txt','.dat')):
        lss_maps_ = open(lss_maps_,'r')
        lss_maps_ = lss_maps_.read().splitlines()
    
    if type(lss_masks_) == str and lss_masks_.endswith(('.txt','.dat')):
        lss_masks_ = open(lss_masks_,'r')
        lss_masks_ = lss_masks_.read().splitlines()

    logger.info('Starting Cross-Correlation analysis...')
    cl_txt_f = os.path.join(X_OUT,'%s_%s_crosscorrelation.txt'%(out_label, bin_label))
    if os.path.exists(cl_txt_f):
        logger.info('ATT: Output file already exists!')
        logger.info(cl_txt_f)
        if not kwargs['overwrite']:
            sys.exit()
        else:
            logger.info('ATT: Overwriting files!')
            
    cl_txt = open(cl_txt_f, 'w')

    for i, m1_f in enumerate(fermi_maps_):
        logger.info('Considering Fermi map: %s'%m1_f)
        fermi_map = hp.read_map(m1_f)
        
        if len(fermi_masks_) == len(fermi_maps_):
            fermi_mask_f = fermi_masks_[i]
        elif len(fermi_masks_) == 1:
            fermi_mask_f = fermi_masks_[0]
        else:
            logger.info('ERROR: Inconsistent number of masks!'+
                         'It should be either 1 or the same number as Fermi maps!')
            sys.exit()
        logger.info('Considering Fermi mask: %s'%fermi_mask_f)
        fermi_mask = hp.read_map(fermi_mask_f)
        
        if len(lss_maps_) == len(fermi_maps_):
            m2_f = lss_maps_[i]
        elif len(lss_maps_) == 1:
            m2_f = lss_maps_[0]
        else:
            logger.info('ERROR: Inconsistent number of LSS maps!'+
                         'It should be either 1 or the same number as Fermi maps!')
            sys.exit()
        logger.info('Considering LSS map: %s'%m2_f)
        lss_map = hp.read_map(m2_f)
            
        if len(lss_masks_) == len(fermi_maps_):    
           lss_mask_f = lss_masks_[i]
        elif len(lss_masks_) == 1:
            lss_mask_f = lss_masks_[0]
        else:
            logger.info('ERROR: Inconsistent number of LSS masks!'+
                         'It should be either 1 or the same number as Fermi maps!')
            sys.exit()
        logger.info('Considering LSS mask: %s'%lss_mask_f)
        lss_mask = hp.read_map(lss_mask_f)
        
        logger.info('Computing the pixel window function ...')
        if not len(fermi_map) == len(lss_map):
            logger.info('ERROR: Fermi and LSS maps have not the same NSIDE!')
            sys.exit()
        nside = hp.npix2nside(len(fermi_map))
        wpix = hp.wpix = hp.sphtfunc.pixwin(nside, lmax=l_max-1)
		
        logger.info('Computing the fermi beam window function ...')
        spec = get_powerlaw_spline(gamma)
        emin_str = re.search(r'%s'%EMIN_STR, m1_f)
        emin = float(emin_str.group(0).replace('-',''))
        emax_str = re.search(r'%s'%EMAX_STR, m1_f)
        emax = float(emax_str.group(0).replace('.fits', ''))
        wb_spline = get_1D_wbeam(fermi_wb_matrix, spec, e_min=emin, e_max=emax)
        wb = wb_spline(np.arange(l_max))
        if kwargs['show']:
            wb_spline.plot(show=False, label='%.2f-%.2f GeV'%(emin/1000, emax/1000), color='0.3')
            plt.ylabel('W$_{beam}$', size=15)
            plt.xlabel('$l$', size=15)
            plt.legend(loc=1, fontsize=14)
            plt.grid('off')
            plt.show()
            
        # ------- TO BE DEFINED ------
        if lss_wb_ is not None:
            logger.info('Getting the LSS beam window function ...')
            lss_wb = np.ones(len(wb))
            
        cl_txt.write('MAP1 ---> %s\n'%(m1_f))
        cl_txt.write('MAP2 ---> %s\n'%(m2_f))
        cl_txt.write('ENERGY\t %.2f %.2f \n'%(emin, emax))

        out_folder =  os.path.join(X_OUT, 'output_pol')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        pol_dict = data.POLCEPICE_DICT
        for key in pol_dict:
            if key == 'clfile':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_cl.txt'%(out_label, i))
            if key == 'corfile':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_ccf.txt'%(out_label, i))
            if key == 'cl_outmap_file':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_clraw.txt'%(out_label, i))
            if key == 'covfileout':
                pol_dict[key] = os.path.join(out_folder,'%s_%i_cov.fits'%(out_label, i))
            if key == 'mapfile':
                pol_dict[key] = m1_f
            if key == 'mapfile2':
                pol_dict[key] = m2_f
            if key == 'maskfile':
                pol_dict[key] = fermi_mask_f
            if key == 'maskfile2':
                pol_dict[key] = lss_mask_f
        config_file_name = os.path.join(out_folder, 'pol_config_%s_%i.txt'%(out_label, i))
        
        if lss_wb_ is None:
            wl = np.sqrt(wb * wpix * wpix)
        else:
            # ------- TO BE DEFINED ------
            wl = np.sqrt(wb * lss_wb * wpix * wpix)
        
        if os.path.exists(os.path.join(out_folder,'%s_%i_cl.txt'%(out_label, i))):
            logger.info('Retrivnig APS and Covariance matrix ...')
            _l, _cl, _cl_err = pol_cl_parse(os.path.join(out_folder,'%s_%i_cl.txt'%(out_label, i)),
                                            os.path.join(out_folder,'%s_%i_cov.fits'%(out_label, i)),
                                            wl_array = wl, rebin=True, lmin=lmin, lmax=lmax,
                                            nbin=bin_num, bin_type=bin_alg,
                                            custom_bins=bin_custom, show=kwargs['show'])
            _cov = pol_cov_parse(os.path.join(out_folder,'%s_%i_cov.fits'%(out_label, i)),
                                              wl_array = wl, rebin=True,
                                              nbin=bin_num, bin_type=bin_alg,
                                              lmin=lmin, lmax=lmax,
                                              custom_bins=bin_custom, 
                                              show=kwargs['show'])                      
        else:
            _l, _cl, _cl_err, _cov = pol_cl_calculation(pol_dict, config_file_name, 
                                                        wl_array=wl, rebin=True,
                                                        nbin=bin_num, bin_type=bin_alg,
                                                        lmin=lmin, lmax=lmax,
                                                        custom_bins=bin_custom, 
                                                        show=kwargs['show'])
        cl_txt.write('multipole\t%s\n'%str(list(_l)).replace('[','').replace(']','').replace(', ', ' '))
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[','').replace(']','').replace(', ', ' '))
        cl_txt.write('Cl_ERR\t%s\n'%str(list(_cl_err)).replace('[','').replace(']','').replace(', ', ' '))
        cl_txt.write('COV_FILE ---> %s\n\n\n'%(os.path.join(out_folder,'%s_%i_cov.pkl'%(out_label, i))))
        
    cl_txt.close()

    logger.info('Created %s'%cl_txt_f)



if __name__ == '__main__':
    args = PARSER.parse_args()
    mkCross(**args.__dict__)
