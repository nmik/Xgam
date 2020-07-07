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
from Xgam import X_OUT

OUT_LABEL = 'Test_UGRB-CMB'

# ----> Fermi LAT products

FERMI_MAPS_LIST = [os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_524-1000.fits'),
				   os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_1000-1737.fits'),
                   os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_1737-2754.fits'),
                   os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_2754-4786.fits'),
                   os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_4786-8317.fits'),
                   os.path.join(X_OUT, 'output_flux/10years_SV_t56_gp25src4fgl_glliemv07_flux_8317-14454.fits')]        
FERMI_MASKS_LIST = [os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits'),
 					os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits'),
					os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits'),
					os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits'),
					os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits'),
					os.path.join(X_OUT,'fits/Mask_gp25_4fgl_extsrc.fits')]
					
FERMI_WBEAM_MATRIX = os.path.join(X_OUT, 'P8R3_SOURCEVETO_V2_evt56_Wbeam.txt')
GAMMA = 2.3 #index of the power law assumed to compute the integral of the Wbeam in the energy bin

# ----> LSS tracer products

LSS_TRACER_MAPS_LIST = ['output/fits/CMB_lens_hp512.fits']
LSS_TRACER_MASK_LIST = ['output/fits/mask_lens_hp512.fits']
LSS_TRACER_WBEAM_LIST = [ ]

# ----> APS parameters

MAX_APS_MULTIPOLE = 1500
BINNING_MIN_MULTIPOLE = 0
BINNING_MAX_MULTIPOLE = 1500
BINNING_MULTIPOLE_NBIN = 21
BINNING_MULTIPOLE_ALGORITHM = 'lin'  #'lin' or 'log'
BINNING_CUSTOM = None #None or list. If a list, previous parameters are neglected.
BINNING_LABEL = '20bins'


# ----> PolSpice dictionary 

POLCEPICE_DICT = {'mapfile' : 'DEFAULT',
                  'mapfile2' : 'DEFAULT',
                  'maskfile' : 'DEFAULT',
                  'maskfile2' : 'DEFAULT',
                  'clfile' : 'DEFAULT',
                  'cl_outmap_file' : 'DEFAULT',
                  'covfileout' : 'DEFAULT',
                  'corfile' : 'DEFAULT',
				  'apodizesigma' : 100,
                  'thetamax' : 100,
                  'nlmax' : MAX_APS_MULTIPOLE,
                  'verbosity' : 1,
                  'apodizetype' : 'NO',
                  'beam' : 'NO',
                  'beam_file' : 'NO',
                  'beam2' : 'NO',
                  'beam_file2' : 'NO',
                  'cl_inmap_file' : 'NO',
                  'cl_outmask_file' : 'NO',
                  'cl_inmask_file' : 'NO',
                  'decouple' : 'NO',
                  'dry' : 'NO',
                  'extramapfile' : 'NO',
                  'extramapfile2' : 'NO',
                  'fits_out' : 'NO',
                  'kernelsfileout' : 'NO',
                  'listmapfiles1_1' : 'NO',
                  'listmapfiles1_2' : 'NO',
                  'listmapfiles1_3' : 'NO',
                  'listmapfiles1_4' : 'NO',
                  'listmapfiles1_5' : 'NO',
                  'listmapfiles1_6' : 'NO',
                  'listmapfiles1_7' : 'NO',
                  'listmapfiles1_8' : 'NO',
                  'listmapfiles1_9' : 'NO',
                  'listmapfiles1_10' : 'NO',
                  'listmapfiles2_1' : 'NO',
                  'listmapfiles2_2' : 'NO',
                  'listmapfiles2_3' : 'NO',
                  'listmapfiles2_4' : 'NO',
                  'listmapfiles2_5' : 'NO',
                  'listmapfiles2_6' : 'NO',
                  'listmapfiles2_7' : 'NO',
                  'listmapfiles2_8' : 'NO',
                  'listmapfiles2_9' : 'NO',
                  'listmapfiles2_10' : 'NO',
                  'listmapweights1_1' :  1.00000000000000,
                  'listmapweights1_2' : 1.00000000000000,
                  'listmapweights1_3' : 1.00000000000000,
                  'listmapweights1_4' : 1.00000000000000,
                  'listmapweights1_5' : 1.00000000000000,
                  'listmapweights1_6' : 1.00000000000000,
                  'listmapweights1_7' : 1.00000000000000,
                  'listmapweights1_8' : 1.00000000000000,
                  'listmapweights1_9' : 1.00000000000000,
                  'listmapweights1_10' : 1.0000000000000,
                  'listmapweights2_1' : 1.00000000000000,
                  'listmapweights2_2' : 1.00000000000000,
                  'listmapweights2_3' : 1.00000000000000,
                  'listmapweights2_4' : 1.00000000000000,
                  'listmapweights2_5' : 1.00000000000000,
                  'listmapweights2_6' : 1.00000000000000,
                  'listmapweights2_7' : 1.00000000000000,
                  'listmapweights2_8' : 1.00000000000000,
                  'listmapweights2_9' : 1.00000000000000,
                  'listmapweights2_10' : 1.00000000000000,
                  'maskfilep' : 'YES',
                  'maskfilep2' : 'YES',
                  'normfac' : 'NO',
                  'npairsthreshold' : 'NO',
                  'noisecorfile' : 'NO',
                  'noiseclfile' : 'NO',
                  'overwrite' : 'YES',
                  'polarization' : 'NO',
                  'pixelfile' : 'NO',
                  'subav' : 'NO',
                  'subdipole' : 'NO',
                  'symmetric_cl' : 'NO',
                  'tf_file' : 'NO',
                  'tolerance' : 'NO',
                  'weightfile' : 'NO',
                  'weightfilep' : 'YES',
                  'weightfile2' : 'NO',
                  'weightfilep2' : 'YES',
                  'weightpower' : 1.00000000000000,
                  'weightpower2' : 1.00000000000000,
                  'weightpowerp' : 1.00000000000000,
                  'weightpowerp2' : 1.00000000000000,
                  'windowfilein' : 'NO',
                  'windowfileout' : 'NO',
                  }
