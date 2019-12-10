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


import os
import numpy as np
import healpy as hp
from astropy.io import fits as pf
from Xgam.utils.logging_ import logger
import matplotlib.pyplot as plt
from Xgam.utils.matplotlib_ import *

def mask_src(CAT_FILE, MASK_S_RAD, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a                                                                                             
       certain radius away from that point.                                                                                                                           
                                                                                                                                                                      
       cat_file: str                                                                                                                                                  
           .fits file of the sorce catalog                                                                                                                            
       MASK_S_RAD: float                                                                                                                                              
           radius around each source definig bad pixels to mask                                                                                                       
       NSIDE: int                                                                                                                                                     
           healpix nside parameter                                                                                                                                    
    """
    logger.info('Mask for sources activated')
    src_cat = pf.open(CAT_FILE)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    RADrad = np.radians(MASK_S_RAD)
    for i in range (0,len(SOURCES)-1):
        GLON = SOURCES.field('GLON')[i]
        GLAT = SOURCES.field('GLAT')[i]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix)
    BAD_PIX_inrad = []
    for bn in BAD_PIX_SRC:
        pixVec = hp.pix2vec(NSIDE,bn)
        radintpix = hp.query_disc(NSIDE, pixVec, RADrad)
        BAD_PIX_inrad.extend(radintpix)
    BAD_PIX_SRC.extend(BAD_PIX_inrad)
    src_cat.close()
    return BAD_PIX_SRC

def mask_extsrc(CAT_FILE, MASK_S_RAD, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a                                                                                             
       certain radius away from that point.                                                                                                                           
                                                                                                                                                                      
       cat_file: str                                                                                                                                                  
           .fits file of the sorce catalog                                                                                                                            
       MASK_S_RAD: float                                                                                                                                              
           radius around each source definig bad pixels to mask                                                                                                       
       NSIDE: int                                                                                                                                                     
           healpix nside parameter                                                                                                                                    
    """
    logger.info('Mask for extended sources activated')
    src_cat = pf.open(CAT_FILE)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT_EXTENDED = src_cat['ExtendedSources']
    BAD_PIX_SRC = []
    EXT_SOURCES = CAT_EXTENDED.data
    src_cat.close()
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES.field('Source_Name')[i]
        GLON = EXT_SOURCES.field('GLON')[i]
        GLAT = EXT_SOURCES.field('GLAT')[i]
        if 'LMC' in NAME or 'CenA Lobes' in NAME:
        	logger.info('Masking %s with 10 deg radius disk...'%NAME)
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10))
        	BAD_PIX_SRC.extend(radintpix)
        else:
        	logger.info('Masking %s with 5 deg radius disk...'%NAME)
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5))
        	BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC
    
def mask_gp(MASK_GP_LAT, NSIDE):
    """Returns the 'bad pixels' around the galactic plain .                                                                                                           
                                                                                                                                                                      
       MASK_GP_LAT: float                                                                                                                                             
           absolute value of galactic latitude definig bad pixels to mask                                                                                             
       NSIDE: int                                                                                                                                                     
           healpix nside parameter                                                                                                                                    
    """
    logger.info('Mask for the galactic plane activated')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_GP = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if abs(b) <= MASK_GP_LAT:
            BAD_PIX_GP.append(iii[i])
    return BAD_PIX_GP

