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


"""Logging utilities, building on top of the python logging module.
"""

import logging
import sys
import time


logger = logging.getLogger('Xgam')
logger.setLevel(logging.INFO)


""" Configure the main terminal logger.
"""
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.INFO)
consoleFormatter = logging.Formatter(">>> %(message)s")
consoleHandler.setFormatter(consoleFormatter)
logger.addHandler(consoleHandler)


def abort(message = ''):
    """Abort the execution (via a sys.exit) with a message.
    Use this with care, and opt for custom exceptions whenever possible.
    """
    if message != '':
        message = '%s. Abort.' % message
    else:
        message = 'Abort.'
    sys.exit(message)

def startmsg():
    """Print the start message.
    """
    BUILD_DATE = '2019-2020'
    FERMI_DATA_RELEASE = 'P8R3_V2'
    print('\n    Welcome to Xgam (built on %s).\n' %BUILD_DATE)
    print('    Autors: Michela Negro, NASA-GSFC/CRESST-UMBC.')
    print('    On behalf of the Fermi-LAT Collaboration.\n')
    print('    This is a framework created to prepare Fermi-LAT data')
    print('    for some cross-correlation analyses. \n')
    print('    (Fermi data release %s). \n\n' %FERMI_DATA_RELEASE)

if __name__ == '__main__':
    startmsg()
