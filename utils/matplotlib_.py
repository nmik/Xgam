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


"""matplotlib configuration module.
"""

import matplotlib
import numpy as np
import time
import os
from matplotlib import pyplot as plt

from Xgam import X_OUT_FIG
from Xgam.utils.logging_ import logger

DEFAULT_FIG_WIDTH = 8.
DEFAULT_FIG_HEIGHT = 6.
DEFAULT_FIG_SIZE = (DEFAULT_FIG_WIDTH, DEFAULT_FIG_HEIGHT)

def property(key):
    """Return a given matplotlib configuration property.
    """
    return matplotlib.rcParams[key]

def context_two_by_two(scale=1.9):
    """Setup the current figure for a 2x2 panel.
    """
    _size = (scale*DEFAULT_FIG_WIDTH, scale*DEFAULT_FIG_HEIGHT)
    _rc = {'figure.figsize': _size}
    return matplotlib.rc_context(rc = _rc)

def context_no_grids():
    """Setup the current figure with no grids.
    """
    _rc = {'axes.grid': False}
    return matplotlib.rc_context(rc = _rc)

def overlay_tag(x=0.95, y=0.95, color='black'):
    """Overlay the Xcorr tag on the current figure.
    """
    text = 'Created by Xgam on %s' % (time.asctime())
    plt.text(x, y, text, color=color, size=10, horizontalalignment='right',
                transform=plt.gca().transAxes)

def setup():
    """Basic setup.
    """
    matplotlib.rc('figure', facecolor='white', figsize=DEFAULT_FIG_SIZE)
    matplotlib.rc('lines', linewidth=2.0)
    matplotlib.rc('patch', linewidth=0.5, facecolor='blue', edgecolor='eeeeee',
                  antialiased=True)
    matplotlib.rc('mathtext', fontset='cm')
    matplotlib.rc('axes',facecolor='white', edgecolor='bcbcbc', grid=False,
                  titlesize=18, labelsize=16)
    matplotlib.rc('xtick', labelsize=14)
    matplotlib.rc('ytick', labelsize=14)
    matplotlib.rc('legend', fancybox=True, frameon=False, numpoints=1, fontsize=16)

setup()

if __name__ == '__main__':

    """Unit test for matplotlib setup
    """
    x = np.logspace(-1, 3, 1000)
    y1 = 100*x**(-2.7)
    y2 = 10*x**(-1.5)
    plt.figure()
    plt.plot(x, y1, label='PL 1')
    plt.plot(x, y2, label='PL 2')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Flux [cm$^{-2}$s$^{-1}$]')
    plt.title('Test plot settings')
    plt.legend()
    overlay_tag()
    plt.show()
