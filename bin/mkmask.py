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

Last login: Tue Dec 10 08:32:12 on ttys003
(base) gs66-bloodymary:~ mnegro$ ssh -XY mnegro@zoroastro.to.infn.it
mnegro@zoroastro.to.infn.it's password: 
Last login: Tue Nov 19 22:31:15 2019 from wcne-128-154-203-56.gsfc.nasa.gov
Linux to05xl.to.infn.it 2.6.32-504.el6.x86_64 #1 SMP Tue Sep 16 01:56:35 EDT 2014 x86_64 x86_64 x86_64 GNU/Linux

===============================================================================

                      T O 0 5 X L . T O . I N F N . I T 
                         (aka Zoroastro/Zarathustra)

                          This is a private system.
           Do not attempt to login unless you are an authorized user.
       Any authorized or unauthorized access and use may be monitored and
       can result in criminal or civil prosecution under applicable law.

                    Machines available for public login:

              Linux: TO05XL [aka ZOROASTRO]       (RHEL6 x86_64)

                     TO01XL [aka ARTABAN]         (RHEL6 x86_64)
                     TO02XL [aka MELCHIORRE]      (RHEL6 x86_64)
                     TO03XL [aka GASPARE]         (RHEL6 x86_64)
                     TO04XL [aka BALDASSARRE]     (RHEL6 x86_64)

              Tru64: TOMMASO   (5.1B  Faster access to /scratch)
                     PATRIZIA  (5.1B-4)

===============================================================================

stty: invalid argument `new'
Try `stty --help' for more information.
to05xl.to.infn.it> ssh -XY mnegro@fermi
mnegro@fermi's password: 
Last login: Tue Nov 19 22:32:13 2019 from to05xl.to.infn.it
[mnegro@fermi ~]$ cd /data1/software/
Anaconda2-2019.03-Linux-x86_64.sh  mnegro_aniso_software/
Anaconda2-4.3.1-Linux-x86_64.sh    
[mnegro@fermi ~]$ cd /data1/software/mnegro_aniso_software/GRATools/
[mnegro@fermi GRATools]$ ll
total 96K
-rw-r--r-- 1 mnegro ater 2.1K Jan  9  2018 __init__.py
drwxr-xr-x 4 mnegro ater  12K Feb  9  2018 Images
drwxr-xr-x 3 mnegro ater 4.0K Jun  6  2018 perAlex
drwxr-xr-x 2 mnegro ater 4.0K Oct 30  2018 sh
drwxrwxr-x 2 mnegro ater 4.0K Jan  7  2019 bin
drwxrwxr-x 7 mnegro ater  44K Feb 12  2019 output
drwxr-xr-x 2 mnegro ater 4.0K Mar  4  2019 __pycache__
drwxrwxr-x 3 mnegro ater 4.0K Mar  8  2019 utils
drwxrwxr-x 6 mnegro ater  12K Apr 16  2019 config
drwxr-xr-x 2 mnegro ater 4.0K Nov  3 16:28 test
[mnegro@fermi GRATools]$ cd utils/
[mnegro@fermi utils]$ ll
total 384K
-rw-r--r-- 1 mnegro ater    0 Jul 21  2016 __init__.py
-rw-r--r-- 1 mnegro ater 2.1K Jul 22  2016 logging_.py
-rw-r--r-- 1 mnegro ater 3.3K Jul 25  2016 matplotlib_.py
-rw-r--r-- 1 mnegro ater 8.8K May 30  2017 ScienceTools_.py
-rw-r--r-- 1 mnegro ater  24K Jun 17  2017 gSpline.py
-rw-r--r-- 1 mnegro ater    8 Aug  1  2017 #gPolSpice.py]#
-rw-r--r-- 1 mnegro ater  142 May 24  2018 __init__.pyc
-rw-r--r-- 1 mnegro ater 1.6K May 24  2018 logging_.pyc
-rw-r--r-- 1 mnegro ater 3.6K May 24  2018 matplotlib_.pyc
-rw-r--r-- 1 mnegro ater  28K May 24  2018 gSpline.pyc
-rw-r--r-- 1 mnegro ater  12K Sep 17  2018 gDrawRef.py
-rw-r--r-- 1 mnegro ater  13K Sep 17  2018 gDrawRef.pyc
-rw-r--r-- 1 mnegro ater  11K Sep 19  2018 gSourceTemplate.py
-rw-r--r-- 1 mnegro ater 7.6K Oct 24  2018 ScienceTools_.pyc
-rw-r--r-- 1 mnegro ater  16K Oct 25  2018 gFTools.py
-rw-r--r-- 1 mnegro ater  16K Oct 25  2018 gFTools.pyc
-rw-r--r-- 1 mnegro ater 8.9K Oct 26  2018 gSourceTemplate.pyc
-rw-r--r-- 1 mnegro ater  14K Nov  8  2018 gWindowFunc.py
-rw-r--r-- 1 mnegro ater  14K Nov 12  2018 gWindowFunc.pyc
-rw-r--r-- 1 mnegro ater  27K Nov 19  2018 gForeground.py
-rw-r--r-- 1 mnegro ater  27K Nov 19  2018 #gForeground.py#
-rw-r--r-- 1 mnegro ater  20K Dec 13  2018 gMasks.py
-rw-r--r-- 1 mnegro ater  17K Dec 13  2018 gMasks.pyc
-rw-r--r-- 1 mnegro ater  23K Dec 13  2018 gForeground.pyc
-rw-r--r-- 1 mnegro ater  14K Dec 13  2018 gPolSpice.py
-rw-r--r-- 1 mnegro ater  14K Dec 13  2018 gPolSpice.pyc
drwxr-xr-x 2 mnegro ater 4.0K Mar  4  2019 __pycache__
-rw-r--r-- 1 mnegro ater  16K Mar  8  2019 #gFTools.py#
[mnegro@fermi utils]$ emacs -nw gMasks.py


























































[mnegro@fermi utils]$ cd ..
[mnegro@fermi GRATools]$ cd bin/
[mnegro@fermi bin]$ emacs -nw mkmask.py 

File Edit Options Buffers Tools Python Help                                                                                                                           
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

"""Produces masks                                                                                                                                                     
"""

import os
import imp
import numpy as np
import healpy as hp
import pyfits as pf


__description__ = 'Computes fluxes'


"""Command-line switches.                                                                                                                                             
"""
import ast
import argparse
from Xgam.utils.logging_ import logger

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='sources mask activated')
PARSER.add_argument('--extsrcmask', type=ast.literal_eval,
                    choices=[True, False], default=False,
                    help='extended sources mask activated')
PARSER.add_argument('--gpmask', type=str, choices=['no','shape', 'flat'],
                    default='no',
                    help='galactic plain mask (only "flat" available now)')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkMask(**kwargs):
    """Routine to produce sky maps (limited edition)                                                                                                                                                               
    """
    logger.info('Starting mask production...')
    get_var_from_file(kwargs['config'])
    bad_pix = []
    nside = data.NSIDE
    out_label = data.OUT_LABEL
    energy = data.ENERGY
    npix = hp.nside2npix(nside)
    mask = np.ones(npix)
    if kwargs['srcmask'] == True:
        from GRATools.utils.gMasks import mask_src
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        bad_pix += mask_src(cat_file, src_mask_rad, nside)
    if kwargs['extsrcmask'] == True:
        from GRATools.utils.gMasks import mask_extsrc
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        bad_pix += mask_extsrc(cat_file, src_mask_rad, nside)
    if kwargs['gpmask'] == 'flat':
        from GRATools.utils.gMasks import mask_gp
        gp_mask_lat = data.GP_MASK_LAT
        bad_pix += mask_gp(gp_mask_lat, nside)
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0

    out_name = os.path.join(GRATOOLS_CONFIG, 'fits/'+out_label+'.fits')
    fsky = 1-(len(np.unique(bad_pix))/float(npix))
    logger.info('f$_{sky}$ = %.3f'%fsky)
    hp.write_map(out_name, mask, coord='G')
    logger.info('Created %s' %out_name)

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkMask(**args.__dict__)
