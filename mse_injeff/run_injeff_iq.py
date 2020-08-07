#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import time
from mse_injeff.mse_injeff_simu import MsePSF

# Run the code for all field position ...
start = time.time()

# Define parameters
fields = ['X_+0.530', 'X_+0.650', 'X_+0.750', 'Y_+0.000',
          'Y_+0.375', 'Y_+0.530', 'Y_+0.650', 'Y_+0.750',
          'Y_-0.375', 'Y_-0.530', 'Y_-0.650', 'Y_-0.750']
zeniths = ['00', '30', '50', '60']
iqs = [0.3, 0.45, 0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.4, 1.55, 1.7, 1.85, 2.0]


def loop_iq(design='segments'):
    # First apply IQ to all field position and zenith distances
    for iq in iqs:
        for zenith in zeniths:
            for field in fields:
                # Initiate PSF
                psf = MsePSF(zenith=zenith, field=field, seeing=iq, design=design)
                print('Init done: ', time.time() - start)
                # Open Zemax
                psf.open_zemax()
                print('Zemax open: ', time.time() - start)
                # Apply IQ
                psf.apply_iq()  # nparray (5, 11, 512, 512)
                print('Convol done: ', time.time() - start)
    return


# Run loop
loop_iq(design='no_segments')


print('hello')