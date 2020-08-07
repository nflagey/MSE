#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from matplotlib import pyplot as plt
import time
from mse_injeff.mse_injeff_simu import MsePSF
from photutils import centroid_com
from photutils.aperture import CircularAperture, aperture_photometry
import pickle

start = time.time()

# Define parameters
fields = ['X_+0.530', 'X_+0.650', 'X_+0.750', 'Y_+0.000',
          'Y_+0.375', 'Y_+0.530', 'Y_+0.650', 'Y_+0.750',
          'Y_-0.375', 'Y_-0.530', 'Y_-0.650', 'Y_-0.750']
zeniths = ['00', '30', '50', '60']
#  Wavelengths
wave = np.asarray([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784])  # in nm
nwave = len(wave)

# Box size for search
boxsize = 21
# Radius range for search
radii = range(25, 100)

# Output
ee80arcsec = {}

# Start loops
for zenith in zeniths:
    ee80arcsec[zenith] = {}
    for field in fields:
        ee80arcsec[zenith][field] = {}
        # Initiate PSF
        psf = MsePSF(zenith=zenith, field=field)
        print('Init done: ', time.time() - start)
        # Open Zemax
        psf.open_zemax()
        speccube = psf.cube[2]  # Select at focus spectral cube
        for w in range(nwave):
            # Select correct slice
            data = speccube[w]
            # Get some info
            norm = np.sum(data)  # Total flux
            x1, y1 = centroid_com(data)  # Center of Mass
            # Define apertures at CoM and with various sizes
            encnrg = np.zeros((boxsize, boxsize, len(radii)))
            for i in range(boxsize):
                for j in range(boxsize):
                    position = [(x1+i-(boxsize-1)/2, y1+j-(boxsize-1)/2)]
                    apertures = [CircularAperture(position, r=r) for r in radii]
                    phot_table = aperture_photometry(data, apertures)
                    encnrg[i, j, :] = [phot_table['aperture_sum_'+str(r)]/norm for r in range(len(radii))]
            # For each radius, find the maximum encircled energy
            max_as_rad = np.max(encnrg, axis=(0, 1))
            # Now identify at which radius we reach 80% ... and convert into EE80 in arcsec
            ee80arcsec[zenith][field][str(wave[w])] = (radii[0] + np.argmin(np.abs(max_as_rad - 0.8))) * 2 * 0.3 / 106.7
            # Plot
            plt.figure()
            plt.imshow(encnrg[:, :, 10])
            plt.colorbar()
            # Save the dictionary (do it at every step to have something to play with while it's running)
            np.save('results/ee80_zemax.npy', ee80arcsec)
            with open('results/ee80_zemax.pkl', 'wb') as f:
                pickle.dump(ee80arcsec, f, pickle.HIGHEST_PROTOCOL)

            print('hello')