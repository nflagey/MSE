#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

"""
Create target lists for cosmology survey and execute Fiber-to-Target Allocator 
"""

# Imports
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
from mse_alloc import MseFibAlloc
import time
# from numpy.random import default_rng

# Initialize RNG
# rng = default_rng()


# 1a. Create targets list(s) for FTA with truncated KiDS catalog
def create_target_list(p_elg=1., p_lbg=1., p_qso=1., nobs=1):
    """
    This function creates target list files for the Fiber-to-Target Allocator.

    :param p_elg: priority score for the ELGs, either a scalar float to give all ELGs the same score, or a string to force random scores
    :param p_lbg: same as p_elg but for the LBGs
    :param p_qso: same as p_elg but for the QSOs
    :param nobs: number of observations required (more than 1 to test dithering or multiple iterations for instance)
    :return: nothing, it just writes files
    """

    # Read full KiDS catalog from Christophe (7,303,813 entries too)
    hdul = fits.open('TARGETS/kids_dr4_171ra190_23.0r24.5.01Apr2020.fits')
    cat = hdul[1].data

    # Get coordinates
    ra = cat['RAJ2000']
    dec = cat['DECJ2000']
    # Get mask
    mask = cat['MASK']
    # Get magnitudes
    umag, gmag, rmag = cat['MAG_GAAP_u'], cat['MAG_GAAP_g'], cat['MAG_GAAP_r']
    imag, zmag, ymag = cat['MAG_GAAP_i'], cat['MAG_GAAP_z'], cat['MAG_GAAP_y']
    jmag, hmag, kmag = cat['MAG_GAAP_J'], cat['MAG_GAAP_H'], cat['MAG_GAAP_Ks']
    # Get magnitudes errors
    umag_err, gmag_err, rmag_err = cat['MAGERR_GAAP_u'], cat['MAGERR_GAAP_g'], cat['MAGERR_GAAP_r']
    imag_err, zmag_err, ymag_err = cat['MAGERR_GAAP_i'], cat['MAGERR_GAAP_Z'], cat['MAGERR_GAAP_Y']
    jmag_err, hmag_err, kmag_err = cat['MAGERR_GAAP_J'], cat['MAGERR_GAAP_H'], cat['MAGERR_GAAP_Ks']
    # Compute colors
    u_g = umag - gmag
    g_r = gmag - rmag
    r_i = rmag - imag
    r_J = rmag - jmag
    r_H = rmag - hmag
    r_Ks = rmag - kmag
    # Galaxy/star selection (0 for galaxies, 1 for stars ... in between most of the time?)
    sgclass = cat['CLASS_STAR']
    # Selection star
    star = (sgclass > 0.985) & (rmag < 23.5) & (rmag_err < 0.2)
    # for bit in [1, 2, 3]: star &= ((mask & 2 ** bit) == 0)

    # Compute coverage
    area = (np.sin(np.min(dec)/180*np.pi) - np.sin(np.max(dec)/180*np.pi)) * (np.min(ra/180*np.pi) - np.max(ra/180*np.pi))
    area *= (180 / np.pi)**2  # square degrees

    # LBG selection
    color_box_lbg = (u_g > 0.0) & (g_r < 1.0) & ((u_g > 2 * g_r + 0.6) | (g_r < 0))
    sel_lbg = color_box_lbg & (rmag > 23.0) & (rmag < 24.5) & (rmag_err < 0.4)
    # ELG selection
    color_box_elg = (u_g < 0.5) & (g_r > 0.0) & (g_r < 0.3) & (r_i < 0.3)
    sel_elg = color_box_elg & (rmag > 23.0) & (rmag < 23.74) & (rmag_err < 0.4)
    # Ly-a QSO selection
    #    Parameters defined by X. Fan
    c1 = 0.95 * u_g + 0.31 * g_r + 0.11 * r_i
    c2 = 0.07 * u_g - 0.49 * g_r + 0.87 * r_i
    c3 = -0.39 * u_g + 0.79 * g_r + 0.47 * r_i
    #    All criteria together
    color_box_qso = (u_g > 0.3) & (c3 < 0.7 - 0.5 * c1) & (c1 > 0.3) & (g_r < 0.6)
    nir_color = (r_J < 1.0) & (r_J > -0.5) & (r_H > -0.3) & (r_Ks > -0.2)
    psf_object = (sgclass > 0.93)
    sel_qso = color_box_qso & (rmag > 19.0) & (rmag < 24.0) & (rmag_err < 0.2) & psf_object & nir_color

    # apply photometric mask
    for bit in [1, 2, 3]:
        sel_lbg &= ((mask & 2**bit) == 0)
        sel_elg &= ((mask & 2**bit) == 0)
        sel_qso &= ((mask & 2**bit) == 0)

    # number of ELG number of LBG
    print(f"{len(ra[sel_lbg])} LBG total, "
          f"{len(ra[sel_lbg]) / area} LBG per sq.deg, "
          f"{len(ra[sel_lbg]) / area * 1.5} LBG per FoV")
    print(f"{len(ra[sel_elg])} ELG total, "
          f"{len(ra[sel_elg]) / area} ELG per sq.deg, "
          f"{len(ra[sel_elg]) / area * 1.5} ELG per FoV")
    print(f"{len(ra[sel_qso])} QSO total, "
          f"{len(ra[sel_qso]) / area} QSO per sq.deg, "
          f"{len(ra[sel_qso]) / area * 1.5} QSO per FoV")

    # make lists for small and large field
    #   (RA, DEC, Nexp, Nrep, umag, gmag, rmag, imag, zmag, Jmag, Hmag, priority, spriority)
    names = ['PID', 'RAJ2000', 'DECJ2000', 'spectro', 'umag', 'gmag', 'rmag', 'imag', 'zmag', 'Jmag', 'Hmag',
             'Nobsreq', 'Nrepeat', 'priority', 'surveypriority', 'Nobsdone']
    for i in range(2):
        if i == 0:
            # make a list within a 2 by 2 sq.degree box
            sel_lbg = sel_lbg & (ra > 178.) & (ra < 180) & (dec < 1) & (dec > -1)
            sel_elg = sel_elg & (ra > 178.) & (ra < 180) & (dec < 1) & (dec > -1)
            sel_qso = sel_qso & (ra > 178.) & (ra < 180) & (dec < 1) & (dec > -1)
            file = 'large'
        else:
            # make a list within a 200 by 200 sq.arcsec box
            sel_lbg = sel_lbg & (ra > 178.9) & (ra < 179) & (dec < .1) & (dec > -.1)
            sel_elg = sel_elg & (ra > 178.9) & (ra < 179) & (dec < .1) & (dec > -.1)
            sel_qso = sel_qso & (ra > 178.9) & (ra < 179) & (dec < .1) & (dec > -.1)
            file = 'small'

        # Deal with random priorities
        if not np.isreal(p_lbg):
            p_lbg_n = np.random.randint(1, 10, len(sel_lbg[sel_lbg])) * 1.
            # rng.choice([0.25, 0.42, 0.69, 0.6, 1., 1.67, 1.44, 2.4, 4], size=len(sel_lbg[sel_lbg]))
        else:
            p_lbg_n = np.ones_like(range(len(sel_lbg[sel_lbg]))) * p_lbg
        if not np.isreal(p_elg):
            p_elg_n = np.random.randint(1, 10, len(sel_elg[sel_elg])) * 1.
        else:
            p_elg_n = np.ones_like(range(len(sel_elg[sel_elg]))) * p_elg
        if not np.isreal(p_qso):
            p_qso_n = np.random.randint(1, 10, len(sel_qso[sel_qso])) * 1.
        else:
            p_qso_n = np.ones_like(range(len(sel_qso[sel_qso]))) * p_qso

        # create Table for LBG
        PID = ['cosmo_lbg' for i in range(len(sel_lbg[sel_lbg]))]
        lr = ['LR' for i in range(len(sel_lbg[sel_lbg]))]
        data_lbg = Table([PID, ra[sel_lbg], dec[sel_lbg], lr,
                          umag[sel_lbg], gmag[sel_lbg], rmag[sel_lbg], imag[sel_lbg], zmag[sel_lbg],
                          jmag[sel_lbg], hmag[sel_lbg],
                          ra[sel_lbg] * 0 + nobs, ra[sel_lbg] * 0 + 1, ra[sel_lbg] * 0 + p_lbg_n, ra[sel_lbg] * 0 + 1,
                          ra[sel_lbg] * 0],
                         names=names)
        # add ELG
        PID = ['cosmo_elg' for i in range(len(sel_elg[sel_elg]))]
        lr = ['LR' for i in range(len(sel_elg[sel_elg]))]
        data_elg = Table([PID, ra[sel_elg], dec[sel_elg], lr,
                          umag[sel_elg], gmag[sel_elg], rmag[sel_elg], imag[sel_elg], zmag[sel_elg],
                          jmag[sel_elg], hmag[sel_elg],
                          ra[sel_elg] * 0 + nobs, ra[sel_elg] * 0 + 1, ra[sel_elg] * 0 + p_elg_n, ra[sel_elg] * 0 + 1,
                          ra[sel_elg] * 0],
                         names=names)
        # add QSO
        PID = ['cosmo_qso' for i in range(len(sel_qso[sel_qso]))]
        lr = ['LR' for i in range(len(sel_qso[sel_qso]))]
        data_qso = Table([PID, ra[sel_qso], dec[sel_qso], lr,
                          umag[sel_qso], gmag[sel_qso], rmag[sel_qso], imag[sel_qso], zmag[sel_qso],
                          jmag[sel_qso], hmag[sel_qso],
                          ra[sel_qso] * 0 + nobs, ra[sel_qso] * 0 + 1, ra[sel_qso] * 0 + p_qso_n, ra[sel_qso] * 0 + 1,
                          ra[sel_qso] * 0],
                         names=names)

        # concatenate tables and save into file
        vstack([data_lbg, data_elg, data_qso]).write('TARGETS/cosmo_targets_' + file + '_new.csv',
                                                     overwrite=True, format='csv')


# 2. Compute number of 5-mn exposures
def run_etc():
    # Assuming all are the same to begin with to follow the cosmology white paper (all exp.times are 1800 seconds)
    return 0


# 3. Fiber-to-target allocation
def run_fta():
    # Execute FTA
    alloc_nodith = MseFibAlloc(file='cosmo_targets_large_new.csv',
                               doplot=False, meth='fixiter', iternum=1, dither=False, allocfrac=90, spectro='LR')
    # Look at results
    res = Table.read('TARGETS/results.csv', format='csv')

    # Fraction of observed targets for various target priority score
    target_score = np.unique(res['priority'] * res['surveypriority'])
    alloc_frac = [np.sum(res['Nobsdone'][res['priority'] * res['surveypriority'] == ts])
                  / np.sum(res['fov'][(res['fov'] == 1) & (res['priority'] * res['surveypriority'] == ts)])
                  for ts in target_score]
    [print(i, j) for i, j in zip(target_score, alloc_frac)]

    # Fraction of observed targets for various program ID
    target_pid = np.unique(res['PID'])
    alloc_frac = [np.sum(res['Nobsdone'][res['PID'] == tp])
                  / np.sum(res['fov'][(res['fov'] == 1) & (res['PID'] == tp)])
                  for tp in target_pid]
    [print(i, j) for i, j in zip(target_pid, alloc_frac)]

    #plt.scatter(target_score, alloc_frac)
    #plt.show()


# Execute code
# --- random priorities
create_target_list(p_elg='random', p_lbg='random', p_qso='random')
# --- fixed priorities
# create_target_list(p_lbg=1, p_elg=1, p_qso=1)
#create_target_list()

t0 = time.time()
run_fta()
print(time.time()-t0)