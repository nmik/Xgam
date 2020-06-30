#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC                                       #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


import os
import sys
import ast
import argparse
import numpy as np
import healpy as hp
import astropy.io.fits as pf
from scipy.interpolate import interp1d

from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.parsing_ import get_energy_from_fits

__description__ = 'Computes fluxes'

"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--foresub', type=ast.literal_eval, choices=[True, False],
                    default='True',
                    help='galactic foreground subtraction activated')
PARSER.add_argument('--overwrite', type=ast.literal_eval, choices=[True, False],
                    default='False',
                    help='overwrite existing maps')
PARSER.add_argument('--nforefit', type=str, choices=['n', 'nlow', 'nhigh'],
                    default='n',
                    help='galactic foreground normalization: N,lower-end, upper-end')


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

def mkRestyle(**kwargs):
    """ Prepare the data sample to the Anisotropy analysis

    """

    logger.info('Starting the restyling...')
    get_var_from_file(kwargs['config'])
    fore_files_ = data.FORE_FILES_LIST
    macro_bins = data.MACRO_BINS
    out_label = data.OUT_LABEL
    mask_label = data.MASK_LABEL
    fore_label = data.FORE_LABEL
    binning_label = data.BINNING_LABEL
    in_labels_list = data.IN_LABELS_LIST
    mask_files_ = data.MASK_FILE
    micro_bin_file = data.MICRO_BINS_FILE
    cguess_file = data.IGRB_FILE
    bincalc = data.BINCALC

    if type(fore_files_) == str and fore_files_.endswith(('.txt','.dat')):
        fore_files_ = open(fore_files_,'r')
        fore_files_ = fore_files_.read().splitlines()
        
    if type(mask_files_) == str and mask_files_.endswith(('.txt','.dat')):
        mask_files_ = open(mask_files_,'r')
        mask_files_ = mask_files_.read().splitlines()

    overwrite = kwargs['overwrite']
    if kwargs['foresub'] == True:
        outfile_name = os.path.join(X_OUT, '%s_%s_%s_%s_datafluxmaps.txt' \
									%(out_label, mask_label, binning_label, fore_label))
    else:
        outfile_name = os.path.join(X_OUT, '%s_%s_%s_noforesub_datafluxmaps.txt' \
									%(out_label, mask_label, binning_label))
    if os.path.exists(outfile_name):
        logger.info('ATT: Output file already exists!')
        logger.info(outfile_name)
        if not overwrite:
            sys.exit()
        else:
            logger.info('ATT: Overwriting files!')
    outfile = open(outfile_name, 'w')
    if kwargs['foresub'] == True:
        outfile.write('#\tE_MIN\tE_MAX\tE_MEAN\tF_MEAN\tFERR_MEAN\tCN\t'+\
        'FSKY\tFORE_N\tFORE_N_errsx\tFORE_N_errdx\tFORE_C\tFORE_C_errsx\t'+\
        'FORE_C_errdx\n')
        fore_N_, fore_N_errsx_, fore_N_errdx_ = [], [], []
        fore_C_, fore_C_errsx_, fore_C_errdx_ = [], [], []
    else:
        outfile.write(
               '# \t E_MIN \t E_MAX \t E_MEAN \t F_MEAN \t FERR_MEAN \t CN \t FSKY\n')

    CN_, FSKY_ = [], []
    E_MIN_, E_MAX_, E_MEAN_ = [], [], []
    FLUX_, FLUX_ERR_ = [], []
    for i, (minb, maxb) in enumerate(macro_bins):
        logger.info('Considering bins from %i to %i...' %(minb, maxb))
        emin, emax, emean = get_energy_from_fits(micro_bin_file, minbinnum=minb, maxbinnum=maxb)
        E_MIN, E_MAX = emin[0], emax[-1]
        E_MEAN = np.sqrt(emax[0]*emin[-1])
        E_MIN_.append(E_MIN)
        E_MAX_.append(E_MAX)
        E_MEAN_.append(E_MEAN)
        logger.info('Emin, Emax, Emean : %.2f, %.2f, %.2f [MeV]'%(E_MIN, E_MAX, E_MEAN))
        logger.info('Summing in time counts and exposures ...')
        micro_bins = np.arange(minb, maxb+1)
        if type(mask_files_) == list:
            mask_file = mask_files_[i]
        else:
            mask_file = mask_files_
        logger.info('Using mask file: %s'%mask_file)
        mask = hp.read_map(mask_file)
        _unmasked = np.where(mask > 0)[0]
        fsky = float(len(_unmasked)*1./len(mask))
        FSKY_.append(fsky)
        logger.info('>>----> FSKY = %e'%fsky)

        time_sum_cnt_, time_sum_exp_ = [], []

        out_cnt_folder = os.path.join(X_OUT, 'output_count')
        if not os.path.exists(out_cnt_folder):
            os.makedirs(out_cnt_folder)
        for b, mb in enumerate(micro_bins):
            micro_cnt_name = os.path.join(out_cnt_folder, '%s_counts_%i.fits'
											%(out_label, mb))
            micro_exp_name = os.path.join(out_cnt_folder, '%s_exposure_%i.fits'
										  %(out_label, mb))
            if os.path.exists(micro_cnt_name) and os.path.exists(micro_exp_name):
                logger.info('Counts and exposure maps ready! Retrieving them...')
                cnt_map = hp.read_map(micro_cnt_name)
                exp_map = hp.read_map(micro_exp_name)
                time_sum_cnt_.append(cnt_map)
                time_sum_exp_.append(exp_map)
            else:
                logger.info('Getting count and exposure files...')
                t_micro_cnt_maps, t_micro_exp_maps = [], []
                for label in in_labels_list:
                    txt_name = os.path.join(X_OUT, '%s_outfiles.txt' %label)
                    txt = open(txt_name,'r')
                    logger.info('Ref: %s'%label)
                    for line in txt:
                        if 'gtbin' in line:
                            cnt_map = (hp.read_map(line.replace('\n', ''), field=mb))
                            t_micro_cnt_maps.append(cnt_map)
                        if 'gtexpcube2' in line:
                            if bincalc == 'CENTER':
                                emap_mean = hp.read_map(line.replace('\n', ''), field=mb)
                            elif bincalc == 'EDGE':
                                emap = hp.read_map(line.replace('\n', ''), field=range(mb, mb+2))
                                emap_mean = np.sqrt(emap[0]*emap[1])
                            else:
                                logger.info('ATT: Invalid bincalc!')
                                sys.exit()
                            t_micro_exp_maps.append(emap_mean)
                    txt.close()
                logger.info('Summing in time and saving micro cnt and exp maps...')
                micro_cnt_map = np.sum(np.array(t_micro_cnt_maps), axis=0)
                hp.write_map(micro_cnt_name, micro_cnt_map, overwrite=overwrite)
                time_sum_cnt_.append(micro_cnt_map)
                micro_exp_map = np.sum(np.array(t_micro_exp_maps), axis=0)
                hp.write_map(micro_exp_name, micro_exp_map, overwrite=overwrite)
                time_sum_exp_.append(micro_exp_map)

        time_sum_cnt_ = np.array(time_sum_cnt_)
        time_sum_exp_ = np.array(time_sum_exp_)

        logger.info('Computing Poisson noise term...')
        npix = len(time_sum_cnt_[0])
        sr = 4*np.pi/npix
        CN_maps = np.array(time_sum_cnt_/time_sum_exp_**2/sr)
        micro_CN = 0
        for b, mb in enumerate(micro_bins):
            micro_CN = micro_CN + np.mean(CN_maps[b][_unmasked])
        CN_.append(micro_CN)
        logger.info('>>----> CN (white noise) term = %e'%micro_CN)

        # now I have finelly gridded (in energy) summed in time fluxes
        time_sum_flux_ = []
        if kwargs['foresub'] == True:
            logger.info('Foreground subtraction activated...')
            from Xgam.utils.foregroundfit_ import fit_foreground_poisson
            from Xgam.utils.foregroundfit_ import get_fore_integral_flux_map

            micro_fore_map_, micro_CN_ = [], []
            micro_fore_N_, micro_fore_N_errsx_, micro_fore_N_errdx_ = [], [], []
            micro_fore_C_, micro_fore_C_errsx_, micro_fore_C_errdx_ = [], [], []
            logger.info('perfom the fit in each micro enrgy bin...')
            for b, mb in enumerate(micro_bins):
                fore_model_map = get_fore_integral_flux_map(fore_files_, emin[b], emax[b])
                micro_fore_map_.append(fore_model_map)
                cguess_data = np.genfromtxt(cguess_file)
                cguess_interp = interp1d(cguess_data[:,0], cguess_data[:,1])
                
                energy = np.sqrt(emin[b]*emax[b])
                delta_en = emax[b]-emin[b]
                c_guess = cguess_interp(energy)
                n, c, n_sx, n_dx, c_sx, c_dx = \
                                   fit_foreground_poisson(fore_model_map,
														  time_sum_cnt_[b],
														  mask_map=mask,
														  exp=time_sum_exp_[b],
														  n_guess=1.,
														  c_guess=c_guess)
                micro_fore_N_.append(n)
                micro_fore_N_errsx_.append(n_sx)
                micro_fore_N_errdx_.append(n_dx)
                micro_fore_C_.append(c)
                micro_fore_C_errsx_.append(c_sx)
                micro_fore_C_errdx_.append(c_dx)

            fore_N_.append(np.mean(micro_fore_N_))
            fore_N_errsx_.append(np.amin(micro_fore_N_errsx_))
            fore_N_errdx_.append(np.amax(micro_fore_N_errdx_))
            fore_C_.append(np.sum(micro_fore_C_))
            fore_C_errsx_.append(np.sum(micro_fore_C_errsx_))
            fore_C_errdx_.append(np.sum(micro_fore_C_errdx_))

            ### compute the flux
            micro_cnt_foresub_map_ = []
            micro_flx_foresub_map_ = []
            micro_flxerr_foresub_map_ = []
            for b, mb in enumerate(micro_bins):
                micro_time_sum_flux = time_sum_cnt_[b]/time_sum_exp_[b]/sr
                if kwargs['nforefit'] == 'nlow':
                    logger.info('Considering N = lower end')
                    fore_norm_flux = micro_fore_N_errsx_[b]*micro_fore_map_[b]
                    fore_norm_count = fore_norm_flux*time_sum_exp_[b]*sr
                    micro_flux_forsesub = micro_time_sum_flux - fore_norm_flux
                    micro_count_forsesub = time_sum_cnt_[b] - fore_norm_count
                    micro_sqrtcount_forsesub = np.sqrt(time_sum_cnt_[b] - fore_norm_count)
                    micro_fluxerr_forsesub = micro_sqrtcount_forsesub/time_sum_exp_[b]/sr
                elif kwargs['nforefit'] == 'nhigh':
                    logger.info('Considering N = upper end')
                    fore_norm_flux = micro_fore_N_errdx_[b]*micro_fore_map_[b]
                    fore_norm_count = fore_norm_flux*time_sum_exp_[b]*sr
                    micro_flux_forsesub = micro_time_sum_flux - fore_norm_flux
                    micro_count_forsesub = time_sum_cnt_[b] - fore_norm_count
                    micro_sqrtcount_forsesub = np.sqrt(time_sum_cnt_[b] - fore_norm_count)
                    micro_fluxerr_forsesub = micro_sqrtcount_forsesub/time_sum_exp_[b]/sr
                else:
                    fore_norm_flux = micro_fore_N_[b]*micro_fore_map_[b]
                    fore_norm_count = (fore_norm_flux*time_sum_exp_[b]*sr).astype(int)
                    micro_flux_forsesub = micro_time_sum_flux - fore_norm_flux
                    micro_count_forsesub = time_sum_cnt_[b]-fore_norm_count
                    micro_count_forsesub[micro_count_forsesub<0] = 0
                    micro_count_forsesub = time_sum_cnt_[b] - fore_norm_count
                    micro_sqrtcount_forsesub = np.sqrt(micro_count_forsesub)
                    micro_fluxerr_forsesub = micro_sqrtcount_forsesub/time_sum_exp_[b]/sr
                micro_cnt_foresub_map_.append(micro_count_forsesub)
                micro_flx_foresub_map_.append(micro_flux_forsesub)
                micro_flxerr_foresub_map_.append(micro_fluxerr_forsesub)
            time_sum_cnt_ = np.array(micro_cnt_foresub_map_)
            time_sum_flux_ = np.array(micro_flx_foresub_map_)
            time_sum_fluxerr_ = np.array(micro_flxerr_foresub_map_)
        else:
            logger.info('Computing the flux for each micro energy bin...')
            time_sum_flux_ = time_sum_cnt_/time_sum_exp_/sr
            time_sum_fluxerr_ = np.sqrt(time_sum_cnt_)/time_sum_exp_/sr
        
        logger.info('Summing in energies micro fluxes ...')
        time_ene_sum_cnt_ = np.sum(time_sum_cnt_, axis=0)
        time_ene_sum_flux_ = np.sum(time_sum_flux_, axis=0)
        time_ene_sum_fluxerr_ = np.sqrt(np.sum(time_sum_fluxerr_**2, axis=0))

        time_ene_sum_flux_masked = hp.ma(time_ene_sum_flux_)
        time_ene_sum_flux_masked.mask = np.logical_not(mask)

        MACRO_MEAN_FLUX = np.average(time_ene_sum_flux_[_unmasked])
        MACRO_MEAN_FLUX_ERR = np.average(time_ene_sum_fluxerr_[_unmasked])
        FLUX_.append(MACRO_MEAN_FLUX)
        FLUX_ERR_.append(MACRO_MEAN_FLUX_ERR)
        logger.info('>>----> Emin, Emax, Emean : %.2f, %.2f, %.2f [MeV]'%(E_MIN, E_MAX, E_MEAN))
        logger.info('>>----> MEAN FLUX = %.2e [cm-2s-1sr-1]' %MACRO_MEAN_FLUX)
        logger.info('>>----> MEAN FLUX ERR / PIX = %.2e [cm-2s-1sr-1]' %MACRO_MEAN_FLUX_ERR)

        logger.info('Saving macro flux maps...')
        out_flx_folder = os.path.join(X_OUT, 'output_flux')
        if kwargs['foresub'] == True:
            macro_cnt_name = os.path.join(out_cnt_folder, '%s_%s_%s_counts_%i-%i.fits'\
									%(out_label, mask_label, fore_label, E_MIN, E_MAX))
            macro_flx_name = os.path.join(out_flx_folder, '%s_%s_%s_flux_%i-%i.fits'\
									%(out_label, mask_label, fore_label, E_MIN, E_MAX))
            macro_flx_msk_name = os.path.join(out_flx_folder, '%s_%s_%s_fluxmasked_%i-%i.fits'\
									%(out_label, mask_label, fore_label, E_MIN, E_MAX))
        else:
            macro_cnt_name = os.path.join(out_cnt_folder, '%s_%s_counts_%i-%i.fits'\
									%(out_label, mask_label, E_MIN, E_MAX))
            macro_flx_name = os.path.join(out_flx_folder, '%s_%s_noforesub_flux_%i-%i.fits'\
									%(out_label, mask_label, E_MIN, E_MAX))
            macro_flx_msk_name = os.path.join(out_flx_folder, '%s_%s_noforesub_fluxmasked_%i-%i.fits'\
									%(out_label, mask_label, E_MIN, E_MAX))
            
        if not os.path.exists(out_flx_folder):
            os.makedirs(out_flx_folder)
        hp.write_map(macro_cnt_name, time_ene_sum_cnt_, overwrite=overwrite)
        hp.write_map(macro_flx_name, time_ene_sum_flux_, overwrite=overwrite)
        hp.write_map(macro_flx_msk_name, time_ene_sum_flux_masked, overwrite=overwrite)

    logger.info('Writing output file...')
    if kwargs['foresub'] == True:
        for i in range(len(FLUX_)):
            outfile.write('%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2e\t%.2f\t%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2e\n' \
						  %(E_MIN_[i], E_MAX_[i], E_MEAN_[i], FLUX_[i], FLUX_ERR_[i],
						    CN_[i], FSKY_[i], fore_N_[i], fore_N_errsx_[i], fore_N_errdx_[i],
						    fore_C_[i], fore_C_errsx_[i], fore_C_errdx_[i]))
    else:
        for i in range(len(FLUX_)):
            outfile.write('%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2e\t%f\n' \
						  %(E_MIN_[i], E_MAX_[i], E_MEAN_[i], FLUX_[i], FLUX_ERR_[i],
						  CN_[i], FSKY_[i]))

    outfile.close()
    logger.info('Created %s'%outfile_name)
    logger.info('Done!')



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkRestyle(**args.__dict__)
