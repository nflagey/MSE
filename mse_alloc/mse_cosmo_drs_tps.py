#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

"""
Figure out if we can find best telescope pointings given a big list of targets
"""

# Imports
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, vstack
import copy
from bokeh.plotting import figure, save
from astropy import wcs
import healpy as hp
from matplotlib import pyplot as plt


# Read full KiDS catalog from Christophe (7,303,813 entries too)
hdul = fits.open('TARGETS/kids_dr4_171ra190_23.0r24.5.01Apr2020.fits')
cat = hdul[1].data
# Get coordinates
ra = cat['RAJ2000']
dec = cat['DECJ2000']

nside = 64  # 128 corresponds to almost 200,000 cells of 0.21 sq.deg. each
npix = hp.nside2npix(nside)
print(npix)
print(hp.nside2pixarea(nside, degrees=True))

pix = hp.ang2pix(nside, ra, dec, lonlat=True)
hist = plt.hist(pix, bins=np.max(pix)-np.min(pix)+1)

m = np.arange(npix)
hp.mollview(m, title="Mollview image RING")
hp.graticule()


#plt.hist(pix, bins=np.max(pix)-np.min(pix)+1)
#plt.show()