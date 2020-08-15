#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table

# - get IQ and AM values
mkam = Table.read('data_fibsize/mkam_dimm_all.csv', format='csv')
mkam_iq_in = mkam['massdimm_seeing']
mkam_am_in = mkam['massdimm_airmass']
# -- remove unwanted values and convert into numpy.arrays
mkam_iq_out = np.asarray(mkam_iq_in[mkam_iq_in != '#NAME?'], dtype=np.float32)
mkam_am_out = np.asarray(mkam_am_in[mkam_iq_in != '#NAME?'], dtype=np.float32)

plt.hist(mkam_iq_out, bins=30, range=(0, 2))
plt.xlabel('IQ (")')
plt.savefig('all_iq.png')

# - get IQ_500
# -- account for airmass to go back to zenith
iq = mkam_iq_out / mkam_am_out ** (3 / 5)
plt.hist(iq, bins=30, range=(0, 2))
plt.xlabel('IQ AM=1.0 (")')
# -- sort values


# -- remove GL
iq = (iq ** (5 / 3) - 0.289 ** (5 / 3)) ** (3 / 5)
# -- compute r0 at 500nm, with IQ in rad
r0 = 0.98 * 500e-9 / (iq / 180 * np.pi / 3600)
# -- convert MKAM-DIMM values into IQ_500 values with L0 = 30m
iq500 = iq * np.sqrt(1 - 2.183 * (r0 / 30) ** 0.356)

plt.hist(iq500, bins=30, range=(0, 2))
plt.xlabel('IQ 500 (")')

print('hello')