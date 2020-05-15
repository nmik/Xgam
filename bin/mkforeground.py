
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

"""Produces foreground maps (work in progress)
"""

import os
import argparse
import numpy as np
import healpy as hp
from astropy.io import fits as pf

from Xgam import X_OUT
from Xgam.utils.logging_ import logger, startmsg
from Xgam.utils.spline_ import xInterpolatedBivariateSplineLinear

__description__ = 'Converter from Cartesian to healpix format'

"""Command-line switches.
"""

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infile', type=str, required=True,
                    help='input fits file of Fermi Foreground model')
PARSER.add_argument('--nsideout', type=int, default=512,
                    help='NSIDE of the output HEALPix maps')
PARSER.add_argument('-s', '--show', type=bool, required=False, default=False,
                    help='input fits file of Fermi Foreground model')

def foreground_map_convert(**kwargs):
    """Viewer interface for healpix maps
    """
    input_file = kwargs['infile']
    nside_out = kwargs['nsideout']
    if not os.path.exists(input_file):
        abort("Map %s not found!"%input_file)
    frmaps = pf.open(input_file)
    maps_slices = frmaps[0].data
    energy = np.array([x[0] for x in frmaps['ENERGIES'].data])
    nside = 2048
    npix = hp.nside2npix(nside)
    iii = np.arange(npix)
    x,y,z = hp.pix2vec(nside, iii)
    lon_hp, lat_hp = hp.rotator.vec2dir(x,y,z,lonlat=True)
    hp_frmap = np.arange(npix, dtype=np.float64)
    lon_fits = np.arange(len(maps_slices[0][0]))
    nresx = 360./len(lon_fits)
    lon_fits_1 = np.flip(lon_fits[1440:]*nresx-180,0)
    lon_fits = np.append(lon_fits_1, np.flip(lon_fits[:1440]*nresx+180,0))
    lat_fits = np.arange(len(maps_slices[0]))
    lat_fits = lat_fits*nresx-90
    fr_e = []
    out_list = []
    i_start = 8
    for i, en in enumerate(energy[i_start:]):
        i = i+i_start
        out_name = os.path.basename(input_file).replace('.fits','_hp%i_%d.fits'
                                                        %(nside_out, en))
        out_path = os.path.join(X_OUT, 'fits', out_name)
        if os.path.exists(out_path):
        	logger.info('Running map conversion for energy %.2f MeV...'%en)
        	logger.info('Retrieving %s'%out_path)
        else:
        	logger.info('Running map conversion for energy %.2f MeV...'%en)
        	frmap = maps_slices[i]
        	fmt = dict(xname='$l$', xunits='deg', yname='$b$', yunits='deg', zname='Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
        	lon, _indexx = np.unique(lon_fits, return_index=True)
        	lat, _indexy = np.unique(lat_fits, return_index=True)
        	frmap = frmap[:, _indexx]
        	frspline = xInterpolatedBivariateSplineLinear(lon, lat, frmap.T, **fmt)
        	for i, pix in enumerate(hp_frmap):
        		hp_frmap[i] = frspline((lon_hp[i]+360)%360, lat_hp[i])
        	fr_e.append(hp_frmap[12426])
        	hp_frmap_out = hp.pixelfunc.ud_grade(hp_frmap, nside_out,  pess=True)
        	hp.write_map(out_path, hp_frmap_out, coord='G')
        	logger.info('Created map %s'%out_path)
        out_list.append(out_path)
    out_name_list = os.path.basename(input_file).replace('.fits','_hp%i_list.txt'%nside_out)
    out_name_list = os.path.join(X_OUT, 'fits', out_name_list)
    logger.info('Writing list of output files: %s'%out_name_list)
    np.savetxt(out_name_list, out_list, fmt='%s')

    if kwargs['show']:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(energy, fr_e, 'o--')
        plt.xlabel('Energy [MeV]')
        plt.xscale('log')
        plt.ylabel('Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
        plt.show()

    frmaps.close()
    logger.info('Done!')

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    foreground_map_convert(**args.__dict__)
