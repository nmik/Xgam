#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, GSFC/CRESST/UMBC    .                                  #
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
from scipy.integrate import quad
from astropy.io import fits as pf
import matplotlib.pyplot as plt

from Xgam.utils.logging_ import logger



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

def mask_extsrc(CAT_FILE, nside=512):
    """Returns the 'bad pixels' defined by the position of a source and a
       certain radius away from that point.

       cat_file: str
           .fits file of the sorce catalog
       mask_rad: float
           radius around each source definig bad pixels to mask
       nside: int
           healpix nside parameter
    """
    logger.info('Mask for extended sources activated')
    src_cat = pf.open(CAT_FILE)
    NPIX = hp.pixelfunc.nside2npix(nside)
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
        	b_pix= hp.pixelfunc.vec2pix(nside, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(nside, (x, y, z), np.radians(10))
        	BAD_PIX_SRC.extend(radintpix)
        else:
        	logger.info('Masking %s with 5 deg radius disk...'%NAME)
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix = hp.pixelfunc.vec2pix(nside, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(nside, (x, y, z), np.radians(5))
        	BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC

def mask_gp(latitude_cut, nside):
    """
    Returns the 'bad pixels' around the galactic plain .

       
    latitude_cut: float
        absolute value of galactic latitude definig bad pixels to mask
    nside: int
        healpix nside parameter (power of 2)
        
    """

    npix = hp.nside2npix(nside)
    
    iii = np.arange(npix)
    x, y, z = hp.pix2vec(nside, iii)
    lon, lat = hp.rotator.vec2dir(x, y, z, lonlat=True)
    
    filter_lat = (abs(lat) < latitude_cut)   
    bad_pix_idx = iii[filter_lat]

    return list(bad_pix_idx)

def mask_src_fluxPSFweighted_1(CAT_FILE, CAT_EXT_FILE, PSF_SPLINE, ENERGY, NSIDE, APODIZE=False):
    """Returns the 'bad pixels' defined by the position of a source and a
       certain radius away from that point. The radii increase with the
       brightness and rescaled by a factor between 1 and 0.3 shaped as the PSF.

       cat_src_file: str
          .fits file with the source catalog
       cat_extsrc_file: str
          .fits file with the extended sources catalog
       ENERGY: float
          Mean energy of the map to be masked
       NSIDE: int
          healpix nside parameter
       APODIZE: bool
       	  if True the apodization of the mask is applied. The fraction of radius to add
       	  to the masked radius for the apodization is k=2.3.
    """
    src_cat = pf.open(CAT_FILE)
    extsrc_cat = pf.open(CAT_EXT_FILE)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    CAT_EXTENDED = extsrc_cat['ExtendedSources']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    EXT_SOURCES = CAT_EXTENDED.data
    FLUX = np.log10(SOURCES.field('Flux1000'))
    src_cat.close()
    extsrc_cat.close()
    psf_en = PSF_SPLINE(ENERGY)
    flux_min, flux_max = min(FLUX), max(FLUX)
    rad_min = 1.5*psf_en
    rad_max = 3*psf_en
    RADdeg = rad_min + FLUX*((rad_max - rad_min)/(flux_max - flux_min)) -\
        flux_min*((rad_max - rad_min)/(flux_max - flux_min))
    RADrad = np.radians(RADdeg)
    logger.info('Masking the extended Sources')
    logger.info('-> 10 deg around CenA and LMC')
    logger.info('-> 5 deg around the remaining')
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES[i][0]
        GLON = EXT_SOURCES.field('GLON')[i]
        GLAT = EXT_SOURCES.field('GLAT')[i]
        if 'LMC' in NAME or 'CenA Lobes' in NAME:
        	logger.info('Masking %s with 10 deg radius disk...'%NAME)
        	rad = 10
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10))
        	BAD_PIX_SRC.extend(radintpix)
        else:
        	logger.info('Masking %s with 5 deg radius disk...'%NAME)
        	rad = 5
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5))
        	BAD_PIX_SRC.extend(radintpix)
    logger.info('Flux-weighted mask for sources activated')
    for i, src in enumerate(SOURCES):
        GLON = SOURCES.field('GLON')[i]
        GLAT = SOURCES.field('GLAT')[i]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix)
        radintpix = hp.query_disc(NSIDE, (x, y, z), RADrad[i])
        BAD_PIX_SRC.extend(radintpix)
    if APODIZE == True:
    	_apd_ring_pix, _apd_ring_val = [], []
    	k = 2.3 # fraction of radius to apodize and add to the radius
    	for i, src in enumerate(SOURCES):
    		apd_rad = k*RADrad[i]
    		GLON = SOURCES.field('GLON')[i]
    		GLAT = SOURCES.field('GLAT')[i]
    		x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
    		b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
    		mask_disk = hp.query_disc(NSIDE, (x, y, z), RADrad[i])
    		apod_disk = hp.query_disc(NSIDE, (x, y, z), apd_rad)
    		apod_ring_pix = np.setxor1d(apod_disk, mask_disk)
    		apod_ring_vec = hp.pixelfunc.pix2vec(NSIDE, apod_ring_pix)
    		apod_ring_dist = hp.rotator.angdist((x,y,z), apod_ring_vec)
    		_apd_ring_pix.append(apod_ring_pix)
    		ang_x = (np.pi/2. * (apod_ring_dist-RADrad[i]))/apd_rad
    		_apd_ring_val.append(np.cos(np.pi/2.-ang_x))
    	return BAD_PIX_SRC, _apd_ring_pix, _apd_ring_val
    else:
        return BAD_PIX_SRC

def compute_flux(E_MIN, E_MAX, SOURCE):
    SPEC_TYPE = SOURCE['SpectrumType']
    E0 = SOURCE['Pivot_Energy']
    if SPEC_TYPE == 'PowerLaw':
        K =  SOURCE['PL_Flux_Density']
        GAMMA = SOURCE['PL_Index']
        dNdE = lambda E: K*(E/E0)**(-GAMMA)

    elif SPEC_TYPE == 'LogParabola':
        K =  SOURCE['LP_Flux_Density']
        alpha = SOURCE['LP_Index']
        beta = SOURCE['LP_beta']
        dNdE = lambda E: K*(E/E0)**(-alpha-beta*np.log(E/E0))

    elif SPEC_TYPE == 'PLSuperExpCutoff':
        K =  SOURCE['PLEC_Flux_Density']
        GAMMA = SOURCE['PLEC_Index']
        a = SOURCE['PLEC_Expfactor']
        b = SOURCE['PLEC_Exp_Index']
        dNdE = lambda E: K*(E/E0)**(-GAMMA)*np.exp(a*(E0**b - E**b))

    return quad(dNdE,E_MIN,E_MAX)[0]

def mask_src_fluxPSFweighted_2(CAT_FILE, CAT_EXT_FILE, PSF_SPLINE, E_MIN, E_MAX, NSIDE, APODIZE=False):
    """Returns the 'bad pixels' defined by the position of a source and a
       certain radius away from that point. The radii increase with the
       brightness and rescaled by a factor between 1 and 0.3 shaped as the PSF.

       cat_src_file: str
            .fits file with the source catalog
       cat_extsrc_file: str
            .fits file with the extended sources catalog
       E_MIN: float
            Lower energy of the map to be masked
       E_MIN: float
            Upper energy of the map to be masked
       NSIDE: int
            healpix nside parameter
       APODIZE: bool
       	    if True the apodization of the mask is applied. The fraction of radius to add
       	    to the masked radius for the apodization is k=2.3.
    """
    src_cat = pf.open(CAT_FILE)
    extsrc_cat = pf.open(CAT_EXT_FILE)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    CAT_EXTENDED = extsrc_cat['ExtendedSources']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    EXT_SOURCES = CAT_EXTENDED.data
    FLUX = []
    for SRC in SOURCES:
        FLUX.append(compute_flux(E_MIN, E_MAX, SRC))
    src_cat.close()
    extsrc_cat.close()
    ENERGY = np.sqrt(E_MIN*E_MAX)
    psf_en = PSF_SPLINE(ENERGY)
    flux_min, flux_max = min(FLUX), max(FLUX)
    FLUX_RATIO = np.array(FLUX)/flux_min
    RADdeg = psf_en*np.sqrt(2*np.log10(5*FLUX_RATIO))
    RADrad = np.radians(RADdeg)
    logger.info('Masking the extended Sources')
    logger.info('-> 10 deg around CenA and LMC')
    logger.info('-> 5 deg around the remaining')
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES[i][0]
        GLON = EXT_SOURCES.field('GLON')[i]
        GLAT = EXT_SOURCES.field('GLAT')[i]
        if 'LMC' in NAME or 'CenA Lobes' in NAME:
        	logger.info('Masking %s with 10 deg radius disk...'%NAME)
        	rad = 10
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10))
        	BAD_PIX_SRC.extend(radintpix)
        else:
        	logger.info('Masking %s with 5 deg radius disk...'%NAME)
        	rad = 5
        	x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        	b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        	BAD_PIX_SRC.append(b_pix)
        	radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5))
        	BAD_PIX_SRC.extend(radintpix)
    logger.info('Flux-weighted mask for sources activated')
    for i, src in enumerate(SOURCES):
        GLON = SOURCES.field('GLON')[i]
        GLAT = SOURCES.field('GLAT')[i]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix)
        radintpix = hp.query_disc(NSIDE, (x, y, z), RADrad[i])
        BAD_PIX_SRC.extend(radintpix)
    if APODIZE == True:
    	_apd_ring_pix, _apd_ring_val = [], []
    	k = 2.3 # fraction of radius to apodize and add to the radius
    	for i, src in enumerate(SOURCES):
    		apd_rad = k*RADrad[i]
    		GLON = SOURCES.field('GLON')[i]
    		GLAT = SOURCES.field('GLAT')[i]
    		x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
    		b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
    		mask_disk = hp.query_disc(NSIDE, (x, y, z), RADrad[i])
    		apod_disk = hp.query_disc(NSIDE, (x, y, z), apd_rad)
    		apod_ring_pix = np.setxor1d(apod_disk, mask_disk)
    		apod_ring_vec = hp.pixelfunc.pix2vec(NSIDE, apod_ring_pix)
    		apod_ring_dist = hp.rotator.angdist((x,y,z), apod_ring_vec)
    		_apd_ring_pix.append(apod_ring_pix)
    		ang_x =  (np.pi/2. * (apod_ring_dist-RADrad[i]))/apd_rad
    		_apd_ring_val.append(np.cos(np.pi/2.-ang_x))
    	return BAD_PIX_SRC, _apd_ring_pix, _apd_ring_val
    else:
        return BAD_PIX_SRC

def mask_south(LAT_LINE, NSIDE):
    """Returns the 'bad pixels' around the galactic plain .

       LAT_LINE: float
           value of celestial declination below which to mask
       NSIDE: int
           healpix nside parameter
    """
    print('Mask for the southern hemishere activated')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if b <= LAT_LINE:
            BAD_PIX.append(iii[i])
    return BAD_PIX

def mask_north(LAT_LINE, NSIDE):
    """Returns the 'bad pixels' around the galactic plain .

       LAT_LINE: float
           value of celestial declination below which to mask
       NSIDE: int
           healpix nside parameter
    """
    print('Mask for the northen hemishere activated')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if b >= LAT_LINE:
            BAD_PIX.append(iii[i])
    return BAD_PIX
