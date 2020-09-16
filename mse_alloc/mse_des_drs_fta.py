#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

"""
Create target lists for DES survey and execute Fiber-to-Target Allocator 
"""

# Imports
from matplotlib import pyplot as plt
import os
import fitsio
from astropy.io import ascii, fits
from astropy.table import Table, vstack
import numpy as np
from mse_alloc import MseFibAlloc


# 1a. Create targets list(s) for FTA with truncated KiDS catalog
def create_target_list(p_stars=1., nobs=1):
    """
    This function creates target list file using DES for the Fiber-to-Target Allocator.

    :param p_stars: priority for stars, either a scalar float to give all ELGs the same score, or a string to force random scores
    :param nobs: number of observations required (more than 1 to test dithering or multiple iterations for instance)
    :return:
    """

    # Read DES extracted catalog
    # (id,ra,dec,mag_auto_g,mag_auto_r,mag_auto_i,mag_auto_z,mag_auto_y,galactic_l,galactic_b,
    # class_star_g,class_star_r,class_star_i,class_star_z,class_star_y)
    cat = Table.read('TARGETS/des_dr1_320ra322_m1dec1.csv')
    # Get coordinates
    ra = cat['ra']
    dec = cat['dec']
    # Get magnitudes
    gmag = cat['mag_auto_g']
    rmag = cat['mag_auto_r']
    imag = cat['mag_auto_i']
    zmag = cat['mag_auto_z']
    ymag = cat['mag_auto_y']
    # Get star/galaxy class
    gstar = cat['class_star_g']
    rstar = cat['class_star_r']
    istar = cat['class_star_i']
    zstar = cat['class_star_z']
    ystar = cat['class_star_y']
    class_star = (gstar + rstar + istar + zstar + ystar) / 5.
    stars = class_star > 0.90
    galaxies = class_star < 0.1

    # Compute coverage
    area = (np.sin(np.min(dec)/180*np.pi) - np.sin(np.max(dec)/180*np.pi)) * (np.min(ra/180*np.pi) - np.max(ra/180*np.pi))
    area *= (180 / np.pi)**2  # square degrees

    # create Table for LBG
    PID = ['DES' for i in range(len(stars[stars]))]
    data = Table([PID, ra[stars], dec[stars],
                  gmag[stars] * 0, gmag[stars], rmag[stars], imag[stars], zmag[stars],
                  zmag[stars] * 0, zmag[stars] * 0,
                  ra[stars] * 0 + nobs, ra[stars] * 0 + 1, ra[stars] * 0 + p_stars, ra[stars] * 0 + 1,
                  ra[stars] * 0],
                 names=['PID', 'RAJ2000', 'DECJ2000',
                        'umag', 'gmag', 'rmag', 'imag', 'zmag', 'Jmag', 'Hmag',
                        'Nobsreq', 'Nrepeat', 'priority', 'surveypriority', 'Nobsdone'])
    data.write('TARGETS/des_targets.csv', overwrite=True, format='csv')

# 2. Compute number of 5-mn exposures
def run_etc():
    # Assuming all are the same to begin with to follow the cosmology white paper (all exp.times are 1800 seconds)
    return 0


# 3. Fiber-to-target allocation
def run_fta():
    # Execute FTA
    alloc_nodith = MseFibAlloc(file='des_targets.csv',
                               doplot=False, meth='fixiter', iternum=1, dither=False, allocfrac=90, spectro='HR')
    # Look at results
    res = Table.read('TARGETS/results.csv', format='csv')

    # Plot fraction of observed targets for various target score
    target_score = np.unique(res['priority'] * res['surveypriority'])
    alloc_frac = [np.sum(res['Nobsdone'][res['priority'] * res['surveypriority'] == ts])
                  / np.sum(res['fov'][(res['fov'] == 1) & (res['priority'] * res['surveypriority'] == ts)])
                  for ts in target_score]
    print(target_score, alloc_frac)

    #plt.scatter(target_score, alloc_frac)
    #plt.show()


# Execute code
create_target_list()

run_fta()