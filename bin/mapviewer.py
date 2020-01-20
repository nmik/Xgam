
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

"""Visualize the healpix maps
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

from Xgam import X_OUT
from Xgam.utils.matplotlib_ import *
from Xgam.utils.logging_ import logger, startmsg

__description__ = 'Visualize the maps'


"""Command-line switches.
"""
formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infile', type=str, required=True,
                    help='input fits file')
PARSER.add_argument('--udgrade', type=int, default=512,
                    help='apply a up or under grade of the nside')
PARSER.add_argument('--norm', type=str, choices=['lin', 'log','hist','optimized'],
                    default='lin', help='specify the scale of the z axis')
PARSER.add_argument('--coord', type=str, choices=['G','E','C'], default='G',  nargs='*',
					help='Define the coord sys. If two agrs are provided then converts'+
					'the coordinates from the firsts to the seconds.')
PARSER.add_argument('--cmap', type=str, default='viridis', help='color scale')
PARSER.add_argument('--title', type=str, default='Map', help='map title')
PARSER.add_argument('--unit', type=str, default='', help='units to display on the bar')
PARSER.add_argument('--smoothing', type=float, default=None,
                    help='Gaussian beam fwhm for the smoothing (in deg)')
PARSER.add_argument('--minval', type=float, default=None, help='Minimum range value')
PARSER.add_argument('--maxval', type=float, default=None, help='Maximum range value')
PARSER.add_argument('--counts', type=ast.literal_eval, choices=[True, False],
                    default=False, help='set to True if the map is a count map')
PARSER.add_argument('-s', '--save', type=ast.literal_eval, choices=[True, False],
                    default=False, help='set to True to save the plot map')


def maps_view(**kwargs):
    """Viewer interface for healpix maps                                                                                                                                                             
    """
    TITLE = kwargs['title']
    NORM = kwargs['norm']
    MAXVAL = kwargs['maxval']
    MINVAL = kwargs['minval']
    CMAP = kwargs['cmap']
    COORD = kwargs['coord']
    if len(COORD) > 1:
    	COORD = [COORD[0], COORD[1]]
    UNIT = kwargs['unit']
    SMOOTH = kwargs['smoothing']
    input_file = kwargs['infile']
    if not os.path.exists(input_file):
        abort("Map %s not found!"%input_file)
    healpix_maps = hp.read_map(input_file)
    if kwargs['smoothing'] is not None:
    	healpix_maps = hp.sphtfunc.smoothing(healpix_maps, fwhm=np.radians(SMOOTH))
    else:
    	pass
    nside_out = kwargs['udgrade']
    if kwargs['counts'] == False:
        healpix_map = hp.pixelfunc.ud_grade(healpix_maps, nside_out, pess=True)
    else:
    	healpix_map = hp.pixelfunc.ud_grade(healpix_maps, nside_out, pess=True, power=-2)
    hp.mollview(healpix_map, title=TITLE, norm=NORM, coord=COORD, 
    					  min=MINVAL, max=MAXVAL, cmap=CMAP, unit=UNIT)
    if kwargs['save'] == True:
    	if not os.path.exists(os.path.join(X_OUT, 'figs')):
    		os.system('mkdir %s' %os.path.join(X_OUT, 'figs'))
    	out_name = os.path.join(X_OUT, 'figs', TITLE.replace(' ', '_'))
    	plt.savefig(out_name)
    plt.show()
	
	
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    maps_view(**args.__dict__)