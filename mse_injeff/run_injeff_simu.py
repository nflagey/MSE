#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
import time
import pickle
from mse_injeff.mse_injeff_simu import MsePSF

# Run the code for all field position ...
start = time.time()

# Define parameters for IQs
fields = ['X_+0.530', 'X_+0.650', 'X_+0.750', 'Y_+0.000',
          'Y_+0.375', 'Y_+0.530', 'Y_+0.650', 'Y_+0.750',
          'Y_-0.375', 'Y_-0.530', 'Y_-0.650', 'Y_-0.750']
zeniths = ['00', '30', '50', '60']
iqs = [0.3, 0.45, 0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.4, 1.55, 1.7, 1.85, 2.0]
spectro_diams = ['HR0.75']
# spectro_diams = ['LR1.00', 'MR1.00', 'HR0.80', 'LR0.90', 'HR0.70']


def loop_ie(design='segments'):
    # Then compute IE
    ie = {}  # result will be a dictionary of dictionaries
    ie_sig = {}

    for iq in iqs:
        ie["iq" + str(iq)] = {}
        ie_sig["iq" + str(iq)] = {}
        for zenith in zeniths:
            ie["iq" + str(iq)]['ZD' + zenith] = {}
            ie_sig["iq" + str(iq)]['ZD' + zenith] = {}
            for field in fields:
                ie["iq" + str(iq)]['ZD' + zenith][field] = {}
                ie_sig["iq" + str(iq)]['ZD' + zenith][field] = {}
                for spectro_diam in spectro_diams:
                    # Load optim_ie (need to do it often as run_psf_optim might be running in parallel)
                    with open('results/'+design+'_injeff_optim_curve_hr075.pkl', 'rb') as f:
                        ie_optim = pickle.load(f)
                    # Initiate PSF
                    psf = MsePSF(zenith=zenith, field=field, seeing=iq, spectro=spectro_diam[0:2],
                                 fibdiam=float(spectro_diam[2:6]), design=design)
                    # Open Zemax (because we need to normalize properly)
                    psf.open_zemax()
                    print('Zemax open: ', time.time() - start)
                    # IE budget
                    psf.ie_budget()
                    print('IE budget done: ', time.time() - start)
                    # Compute IE (and pass optim IE for plot)
                    psf.compute_ie(ie_optim["iq" + str(iq)]['ZD' + zenith][field][spectro_diam])
                    print('IE curve done: ', time.time() - start)
                    # Save IE curve
                    ie["iq" + str(iq)]['ZD' + zenith][field][spectro_diam] = psf.ie_simu
                    ie_sig["iq" + str(iq)]['ZD' + zenith][field][spectro_diam] = psf.ie_simu_stddev

                # Save the IE dictionary (do it at every field to have something to play with while it's running)
                with open('results/'+design+'_injeff_curve_hr075.pkl', 'wb') as f:
                    pickle.dump([ie, ie_sig], f, pickle.HIGHEST_PROTOCOL)
    return


# Run loop
loop_ie(design='no_segments')

print('hello')