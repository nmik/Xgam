#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


"""Remake of some of the Fermi Science Tools
"""


import os
import time
import gt_apps as my_apps
from Xgam.utils.logging_ import logger
from Xgam import X_OUT
from Xgam import FT_DATA_FOLDER


FT_DATA_OUT = os.path.join(FT_DATA_FOLDER,'output')
if not os.path.exists(FT_DATA_OUT):
    logger.info('Creating output folder of selected data (%s)'%FT_DATA_OUT)
    os.makedirs(os.path.join(FT_DATA_OUT))


def gtselect(label, filter_dict):
    """gtselect from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtselect...')
    LABEL = label
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtselect')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    OUTFILE = os.path.join(OUTPATH, LABEL + '_outofselect.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    else:
        for key in filter_dict:
            my_apps.filter[key] = filter_dict[key]
        my_apps.filter['outfile'] = OUTFILE
        my_apps.filter.run()
        logger.info('Created %s'%OUTFILE)
        logger.info('gtselect --> CPU time spent: %.2f'%time.clock())
        return OUTFILE

def gtmktime(label, maketime_dict):
    """gtmktime from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtmktime...')
    LABEL = label
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtmktime')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    OUTFILE = os.path.join(OUTPATH, LABEL + '_outofmktime.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    else:
        for key in maketime_dict:
            if key == 'outfile':
                if maketime_dict[key] == 'DEFAULT':
                    my_apps.maketime['outfile'] = OUTFILE
                else:
                    my_apps.maketime[key] = maketime_dict[key]
                continue
            my_apps.maketime[key] = maketime_dict[key]
        my_apps.maketime.run()
        logger.info('Created %s'%OUTFILE)
        logger.info('gtmktime --> CPU time spent: %.2f'%time.clock())
        return OUTFILE

def gtbin(label, evtbin_dict):
    """gtbin from Science Tools. 

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtbin...')
    LABEL = label
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtbin')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    OUTFILE = os.path.join(OUTPATH, LABEL + '_outofbin.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    else:
        for key in evtbin_dict:
            if key == 'outfile':
                if evtbin_dict[key] == 'DEFAULT':
                    my_apps.evtbin['outfile'] = OUTFILE
                else:
                    my_apps.evtbin[key] = evtbin_dict[key]
                continue
            my_apps.evtbin[key] = evtbin_dict[key]
        my_apps.evtbin.run()
        logger.info('Created %s'%OUTFILE)
        logger.info('gtbin --> CPU time spent: %.2f'%time.clock())
        return OUTFILE

def gtltcube(label, expcube_dict):
    """gtltcube from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtltcube...')
    LABEL = label
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtltcube')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    OUTFILE = os.path.join(OUTPATH, LABEL + '_outofltcube.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    else:
        for key in expcube_dict:
            if key == 'outfile':
                if expcube_dict[key] == 'DEFAULT':
                    my_apps.expCube['outfile'] = OUTFILE
                else:
                    my_apps.expCube[key] = expcube_dict[key]
                continue
            my_apps.expCube[key] = expcube_dict[key]
        my_apps.expCube.run()
        logger.info('Created %s'%OUTFILE)
        logger.info('gtltcube --> CPU time spent: %.2f'%time.clock())
        return OUTFILE

def gtexpcube2(label, expcube2_dict):
    """gtexpcube2 from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtexpcube2...')
    LABEL = label
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtexpcube2')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    OUTFILE = os.path.join(OUTPATH, LABEL + '_expcube.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    else:
        for key in expcube2_dict:
            if key == 'outfile':
                if expcube2_dict[key] == 'DEFAULT':
                    my_apps.gtexpcube2['outfile'] = OUTFILE
                else:
                    my_apps.gtexpcube2[key] = expcube2_dict[key]
                continue
            my_apps.gtexpcube2[key] = expcube2_dict[key]
        my_apps.gtexpcube2.run()
        logger.info('Created %s'%OUTFILE)
        logger.info('gtexpcube2 --> CPU time spent: %.2f'%time.clock())
        return OUTFILE

def gtpsf(gtpsf_dict):
    """gtpsf from Science Tools

       gtpsf_dict: python dict
          To define all the parameters
    """
    expcube = gtpsf_dict['expcube']
    outfile = gtpsf_dict['outfile']
    irfs = gtpsf_dict['irfs']
    evtype = gtpsf_dict['evtype']
    ra = gtpsf_dict['ra']
    dec = gtpsf_dict['dec']
    emin = gtpsf_dict['emin']
    emax = gtpsf_dict['emax']
    nenergies = gtpsf_dict['nenergies']
    thetamax = gtpsf_dict['thetamax']
    ntheta = gtpsf_dict['ntheta']
    os.system('gtpsf expcube=%s outfile=%s irfs=%s evtype=%i ra=%f dec=%f emin=%e emax=%e nenergies=%i thetamax=%i ntheta=%i' \
                  %(expcube, outfile, irfs, evtype, ra, dec, emin, emax, \
                        nenergies, thetamax, ntheta))

def gtEbindef(ebinning_array, file_name='ebinning.txt'):
    """Produces a fits file defining the enrgy binning to fed gtbin.

           ebinning_array: numpy array
               array in which the energy binnin is defined.
    """
    if not os.path.exists(GD_OUT):
        os.makedirs(GD_OUT)
    txt_file_name = os.path.join(GD_OUT, file_name)
    txt_file = open(txt_file_name, 'w')
    fits_file_name = os.path.join(GD_OUT, 
                                  file_name.replace('.txt', '.fits'))
    for emin, emax in zip(ebinning_array[:-1], ebinning_array[1:]):
        txt_file.write('%.4f %.4f\n'%(emin, emax))
    txt_file.close()
    os.system('gtbindef bintype=E binfile=%s outfile=%s energyunits=MeV' \
                  %(txt_file_name, fits_file_name))
    logger.info('Created %s...'%fits_file_name)
    return fits_file_name

def mergeft1(path_to_files, out_file_name, N1week, Nnweek):
    """creates a .txt file with the list of the FT files to merge.

       path_to_files: str
           path where datat files are stored
       out_file_name: str
           name of the txt output file (created in the same folder of data)
       N1week: int
           number of the starting week
       Nnweek: int
           number of the ending week
    """
    if N1week < 9:
        abort('Invalid number of weeks: the minimun must be > or = to 9')
    if Nnweek > 561:
        abort('Invalid number of weeks: the maximum must be < or = to 561')
    outtxtfile = os.path.join(path_to_files, out_file_name)
    if not os.path.exists(outtxtfile):
        out_file = open(outtxtfile, 'w')
        for i in range(N1week, Nnweek+1):
            if i == 9:
                out_file.write("%s/lat_photon_weekly_w00%i_p305_v001.fits \n" \
                                   %(path_to_files,i))
            if i >= 10 and i <= 99:
                out_file.write("%s/lat_photon_weekly_w0%i_p305_v001.fits \n" \
                                   %(path_to_files,i))
            if i > 99:
                if i == 512:
                    pass
                else:
                    out_file.write("%s/lat_photon_weekly_w%i_p305_v001.fits \n" \
                                       %(path_to_files,i))
        out_file.close()
    logger.info('Created %s...' %outtxtfile)
    return '@'+outtxtfile

def mergeft2(path_to_files, out_file_name, N1week, Nnweek):
    """creates a .txt file with the list of the FT files to merge.

       path_to_files: str
           path where datat files are stored
       out_file_name: str
           name of the txt output file (created in the same folder of data)
       N1week: int
           number of the starting week
       Nnweek: int
           number of the ending week
    """
    if N1week < 9:
        abort('Invalid number of weeks: the minimun must be > or = to 9')
    if Nnweek > 561:
        abort('Invalid number of weeks: the maximum must be < or = to 561')
    outtxtfile = os.path.join(path_to_files, out_file_name)
    if not os.path.exists(outtxtfile):
        out_file = open(outtxtfile, 'w')
        for i in range(N1week, Nnweek+1):
            if i == 9:
                out_file.write("%s/lat_spacecraft_weekly_w00"+
                               "%i_p302_v001.fits \n" %(path_to_files,i))
            if i >= 10 and i <= 99:
                out_file.write("%s/lat_spacecraft_weekly_w0"+
                               "%i_p302_v001.fits \n" %(path_to_files,i))
            if i > 99:
                out_file.write("%s/lat_spacecraft_weekly_w"+
                               "%i_p302_v001.fits \n" %(path_to_files,i))
        out_file.close()
    logger.info('Created %s...' %outtxtfile)
    return '@'+outtxtfile

def main():
    """Test section.
    """
    import numpy as np
    from Xgam import FT_DATA_FOLDER

    PH = 'photon'
    SC = 'spacecraft'
    PH_FOLDER = os.path.join(FT_DATA_FOLDER, PH)
    SC_FOLDER = os.path.join(FT_DATA_FOLDER, SC)
    FT1_FILE = mergeft1(PH_FOLDER, 'FT1_w9-12.txt', 9, 12)
    FT2_FILE = mergeft2(PH_FOLDER, 'FT2_w9-12.txt', 9, 12)
    FILTER_CUT='DATA_QUAL==1&&LAT_CONFIG==1&&LAT_MODE==5&&IN_SAA!=T'+\
               '&&((ABS(ROCK_ANGLE)<52))'
    OUT_FILE_LABEL = 'test'
    EBINNING_FILE = gtbindef(np.logspace(3, 3.6, 11))
    data_filter = {'infile': FT1_FILE,
                   'emin': 1000,
                   'emax': 5000,
                   'zmax': 90,
                   'evclass': 1024,
                   'evtype': 32,
                   'clobber': 'yes'}
    out_gtselect = gtselect(OUT_FILE_LABEL, data_filter)
    data_mktime = {'evfile': out_gtselect,
                   'scfile': FT2_FILE,
                   'filter': FILTER_CUT,
                   'roicut': 'no',
                   'clobber': 'yes'}
    out_gtmktime = gtmktime(OUT_FILE_LABEL, data_mktime)
    data_bin = {'evfile': out_gtmktime,
                'algorithm': 'HEALPIX',
                'scfile': FT2_FILE,
                'hpx_ordering_scheme': 'RING',
                'hpx_order': 7,
                'coordsys': 'GAL',                                  
                'hpx_ebin': 'yes',
                'ebinalg': 'FILE',
                'ebinfile': EBINNING_FILE,
                'clobber': 'yes'}
    out_gtbin = gtbin(OUT_FILE_LABEL, data_bin)


if __name__ == '__main__':
    main()
