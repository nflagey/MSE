#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
from mse_itc.mse_fibsize import *
import numpy as np

# Spectrograph and fiber sizes to consider
fibers = ['LR0.90', 'LR0.95', 'LR1.00']
# Source type
src_type = 'extended'  # 'point' or 'extended'

# Other Parameters
# Load one set of IE curves to identify all IQs and ZDs
ie = np.load(here_path + '/THROUGHPUT/no_segments_injeff_curve_for_itc_' + fibers[0] + '_' + src_type + '.npy',
             allow_pickle=True).flat[0]
# IQ/Seeing values (FWHM, in arcsecs, at 550 nm)
iq_vals = list(ie.keys())
# ZD/airmass values
zd_vals = list(ie[iq_vals[0]].keys())
zd_vals = zd_vals[0:3]

# Enter below the code you want to execute
# comp_snr(iq_vals=iq_vals, zd_vals=zd_vals, fibers=fibers, src_type=src_type)
# plot_snr_ratio_all(iq_vals=iq_vals, zd_vals=zd_vals, fibers=fibers, src_type=src_type)
plot_snr_ratio_iq(iq_vals=iq_vals, zd_vals=zd_vals, fibers=fibers, src_type=src_type)

print("Done!")
