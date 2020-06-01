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


"""Science Tools Analysis App.

   This app runs a chain of Fermi Science Tools:
      1) gtselect
      2) gtmktime
      3) gtbin
      4) gtltcube
      5) gtexpcube2
   If an output file already exists it won't be overwritten, because the
   command won't be run. 
   A configuration file is needed, in which all the ScienceTools parameters 
   you want to set must be declared. See 1yr_st_aniso_config.py for e.g.
"""

import os
import sys
import ast
import argparse
import numpy as np

from Xgam import FT_DATA_FOLDER
from Xgam.utils.ScienceTools_ import mergeft1, mergeft2
from Xgam.utils.logging_ import logger, startmsg


__description__ = 'Run a chain of Science Tools'


"""Command-line switches.
"""


formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--gtltcube', type=ast.literal_eval, choices=[True, False], 
                    default=True,
                    help='False if gtltcube command must not be run')

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
    
def mkSTanalysis(**kwargs):
    """Science Tools analysis chain
    """
    assert(kwargs['config'].endswith('.py'))
    get_var_from_file(kwargs['config'])
    logger.info('Starting ST analysis...')

    from Xgam.utils.ScienceTools_ import mergeft1
    PH, SC = 'photon', 'spacecraft'
    PH_FOLDER = os.path.join(FT_DATA_FOLDER, PH)
    SC_FOLDER = os.path.join(FT_DATA_FOLDER, SC)
    start_week, end_week = data.START_WEEK, data.STOP_WEEK
    logger.info('Taking data from week %i to week %i'%(start_week, end_week))
    FT1_FILE = mergeft1(PH_FOLDER, 'FT1_w%i-%i.txt'%(start_week, end_week), \
                           start_week, end_week)
    #FT2_FILE = mergeft2(SC_FOLDER, 'FT2_w%i-%i.txt'%(start_week, end_week), \
    #                       start_week, end_week)
    out_label = data.OUT_LABEL
    txt_out_files = open('output/'+out_label+'_outfiles.txt', 'w')

    from Xgam.utils.ScienceTools_ import gtselect
    gtselect_dict = data.GTSELECT_DICT
    if gtselect_dict['infile'] == 'DEFAULT':
        gtselect_dict['infile'] = FT1_FILE
    out_gtselect = gtselect(out_label, gtselect_dict)
    txt_out_files.write(out_gtselect+'\n')

    from Xgam.utils.ScienceTools_ import gtmktime
    gtmktime_dict = data.GTMKTIME_DICT
    if gtmktime_dict['evfile'] == 'DEFAULT':
        gtmktime_dict['evfile'] = out_gtselect
    if gtmktime_dict['scfile'] == 'DEFAULT':
        gtmktime_dict['scfile'] = FT2_FILE
    out_gtmktime = gtmktime(out_label, gtmktime_dict)
    txt_out_files.write(out_gtmktime+'\n')

    from Xgam.utils.ScienceTools_ import gtbin
    gtbin_dict = data.GTBIN_DICT
    if gtbin_dict['evfile'] == 'DEFAULT':
        gtbin_dict['evfile'] = out_gtmktime
    if gtbin_dict['scfile'] == 'DEFAULT':
        gtbin_dict['scfile'] = FT2_FILE
    out_gtbin = gtbin(out_label, gtbin_dict)
    txt_out_files.write(out_gtbin+'\n')

    if kwargs['gtltcube'] == True:
        from Xgam.utils.ScienceTools_ import gtltcube
        gtltcube_dict = data.GTLTCUBE_DICT
        if gtltcube_dict['evfile'] == 'DEFAULT':
            gtltcube_dict['evfile'] = out_gtmktime
        if gtbin_dict['scfile'] == 'DEFAULT':
            gtbin_dict['scfile'] = FT2_FILE
        out_gtltcube = gtltcube(out_label, gtltcube_dict)
        txt_out_files.write(out_gtltcube+'\n')
    else:
        logger.info('Not running gtltcube.')
        pass

    from Xgam.utils.ScienceTools_ import gtexpcube2
    gtexpcube2_dict = data.GTEXPCUBE2_DICT
    if gtexpcube2_dict['infile'] == 'DEFAULT':
        gtexpcube2_dict['infile'] = out_gtltcube
    if gtexpcube2_dict['cmap'] == 'DEFAULT':
        gtexpcube2_dict['cmap'] = out_gtbin
    out_gtexpcube2 = gtexpcube2(out_label, gtexpcube2_dict)
    txt_out_files.write(out_gtexpcube2+'\n')
    txt_out_files.close()
    logger.info('Created output/'+out_label+'_outfiles.txt')
    logger.info('Done!')


if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkSTanalysis(**args.__dict__)