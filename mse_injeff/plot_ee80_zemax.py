#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
from matplotlib import pyplot as plt
import numpy as np
import pickle
from matplotlib import colors
from astropy.table import Table
from astropy.io import fits
import copy

# Define parameters
fields = ['Y_+0.000', 'Y_-0.375', 'Y_+0.375', 'Y_-0.530', 'X_+0.530',
          'Y_+0.530', 'Y_+0.650', 'Y_-0.650', 'X_+0.650',
          'Y_+0.750', 'Y_-0.750', 'X_+0.750']
yfields = ['Y_-0.750', 'Y_-0.650', 'Y_-0.530', 'Y_-0.375', 'Y_+0.000', 'Y_+0.375', 'Y_+0.530', 'Y_+0.650', 'Y_+0.750']
xfields = ['X_+0.750', 'X_+0.650', 'X_+0.530', 'Y_+0.000', 'X_+0.530', 'X_+0.650', 'X_+0.750']
ypos = [-.75, -.65, -.53, -.375, 0, 0.375, 0.53, 0.65, 0.75]
xpos = [-.75, -.65, -.53, 0, 0.53, 0.65, 0.75]
zeniths = ['00', '30', '50', '60']
#  Wavelengths
waves = np.asarray([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784])  # in nm
nwave = len(waves)
# Colors
waveHSV = np.zeros([11, 3])
waveHSV[:, 0] = np.arange(nwave, 0, -1) / nwave
waveHSV[:, 1] = 1
waveHSV[:, 2] = 1
wave_color = [colors.hsv_to_rgb(waveHSV[w, :]) for w in range(nwave)]

# Open results
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/ee80/ee80_zemax.pkl', 'rb') as f:
    ee80arcsec = pickle.load(f)
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/ee80/ee80_zemax_monolithic.pkl', 'rb') as f:
    ee80arcsec_mono = pickle.load(f)

# Get EE80 directly from Zemax - monolithic M1
# Read Derrick's files
eezmx = {}
for f in range(len(fields)):
    eezmx[fields[f]] = {}
    for w in range(nwave):
        eezmx[fields[f]][str(waves[w])] = Table.read('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/ee80/EE_lam_'+str(w+1)+'.txt',
                                                     format='ascii.basic', guess=False,
                                                     data_start=14+f*88, data_end=99+f*88, header_start=13)
# Get EE80 for each field position and wavelength
ee80zmx = {}
for field in fields:
    ee80zmx[field] = {}
    for w in range(nwave):
        ee80zmx[field][str(waves[w])] = np.interp(0.8, eezmx[field][str(waves[w])]['Fraction'],
                                                  eezmx[field][str(waves[w])]['Radial_distance']) * 2 / 106.7
# Compare with my values
plt.figure()
plt.plot([ee80zmx[field]['551'] for field in fields])
plt.plot([ee80arcsec_mono['30'][field]['551'] for field in fields])
# Save the dictionary
np.save('results/ee80/ee80zemax.npy', ee80zmx)
with open('results/ee80/ee80zemax.pkl', 'wb') as f:
    pickle.dump(ee80zmx, f, pickle.HIGHEST_PROTOCOL)
# Plot those EE80
for zenith in ['30']:
    plt.figure(figsize=(8, 4.5), dpi=120)
    for w in range(nwave):
        plt.plot(xpos, [ee80zmx[xfield][str(waves[w])] for xfield in xfields],
                 color=wave_color[w], label=waves[w])
        plt.plot(ypos, [ee80zmx[yfield][str(waves[w])] for yfield in yfields], color=wave_color[w])
    plt.legend()
    plt.xlabel('Field position (degree)')
    plt.ylabel('EE80 (arcsec)')
    plt.xlim([-.8, .8])
    plt.ylim([0.1, 0.75])
    plt.title('ZD' + zenith)
    plt.savefig('results/ee80/ee80zemax_Z' + zenith + '.png', bbox_inches='tight')
    plt.close()

# Plot each wavelength as function of field position
for zenith in zeniths:
    plt.figure(figsize=(8, 4.5), dpi=120)
    for w in range(nwave):
        plt.plot(xpos, [ee80arcsec[zenith][xfield][str(waves[w])] for xfield in xfields],
                 color=wave_color[w], label=waves[w])
        plt.plot(ypos, [ee80arcsec[zenith][yfield][str(waves[w])] for yfield in yfields], color=wave_color[w])
    plt.legend()
    plt.xlabel('Field position (degree)')
    plt.ylabel('EE80 (arcsec)')
    plt.xlim([-.8, .8])
    plt.ylim([0.1, 0.75])
    plt.title('ZD' + zenith)
    plt.savefig('results/ee80/ee80_Z' + zenith + '.png', bbox_inches='tight')
    plt.close()

# Plot each wavelength as function of field position (monolithic version)
for zenith in zeniths:
    plt.figure(figsize=(8, 4.5), dpi=120)
    for w in range(nwave):
        plt.plot(xpos, [ee80arcsec_mono[zenith][xfield][str(waves[w])] for xfield in xfields],
                 color=wave_color[w], label=waves[w])
        plt.plot(ypos, [ee80arcsec_mono[zenith][yfield][str(waves[w])] for yfield in yfields], color=wave_color[w])
    plt.legend()
    plt.xlabel('Field position (degree)')
    plt.ylabel('EE80 monolithic (arcsec)')
    plt.xlim([-.8, .8])
    plt.ylim([0.1, 0.75])
    plt.title('ZD' + zenith)
    plt.savefig('results/ee80/ee80_mono_Z' + zenith + '.png', bbox_inches='tight')
    plt.close()


# Compare monolithic and segmented values
# plot vs
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    plt.plot([ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths],
             [ee80arcsec_mono[zenith][field][str(waves[w])] for field in fields for zenith in zeniths],
             '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('EE80 monolithic (arcsec)')
plt.xlim([.1, .9])
plt.ylim([.1, .9])
plt.savefig('results/ee80/ee80_vs.png', bbox_inches='tight')
plt.close()
# plot ratio
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    xx = [ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths]
    yy = [ee80arcsec[zenith][field][str(waves[w])] / ee80arcsec_mono[zenith][field][str(waves[w])]
          for field in fields for zenith in zeniths]
    plt.plot(xx, yy, '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('EE80 segmented / EE80 monolithic')
plt.xlim([.1, .9])
plt.ylim([1, 5])
plt.savefig('results/ee80/ee80_ratio.png', bbox_inches='tight')
plt.close()
# plot difference
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    xx = [ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths]
    yy = [ee80arcsec[zenith][field][str(waves[w])] - ee80arcsec_mono[zenith][field][str(waves[w])]
          for field in fields for zenith in zeniths]
    plt.plot(xx, yy, '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('EE80 segmented - EE80 monolithic')
plt.xlim([.1, .9])
plt.ylim([0, .5])
plt.savefig('results/ee80/ee80_diff.png', bbox_inches='tight')
plt.close()
# plot difference scaled with wavelength
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    xx = [ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths]
    yy = [(ee80arcsec[zenith][field][str(waves[w])] - ee80arcsec_mono[zenith][field][str(waves[w])]) / waves[w]
          for field in fields for zenith in zeniths]
    print(np.mean(yy))
    plt.plot(xx, yy, '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('(EE80 segmented - EE80 monolithic) / wavelength')
plt.xlim([.1, .9])
plt.ylim([0, 1e-3])
plt.savefig('results/ee80/ee80_diff_wavscaled.png', bbox_inches='tight')
plt.close()
# plot quadratic diff
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    xx = [ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths]
    yy = [np.sqrt(ee80arcsec[zenith][field][str(waves[w])]**2 - ee80arcsec_mono[zenith][field][str(waves[w])]**2)
          for field in fields for zenith in zeniths]
    plt.plot(xx, yy, '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('Quad. diff. (EE80 segmented - EE80 monolithic)')
plt.xlim([.1, .9])
plt.ylim([.1, .9])
plt.savefig('results/ee80/ee80_quaddiff.png', bbox_inches='tight')
plt.close()
# plot quadratic diff scaled with wavelength
plt.figure()
plt.plot([0, 1], [0, 1], 'k')
for w in range(nwave):
    xx = [ee80arcsec[zenith][field][str(waves[w])] for field in fields for zenith in zeniths]
    yy = [np.sqrt(ee80arcsec[zenith][field][str(waves[w])]**2 - ee80arcsec_mono[zenith][field][str(waves[w])]**2) / waves[w]
          for field in fields for zenith in zeniths]
    plt.plot(xx, yy, '+', color=wave_color[w], label=waves[w])
plt.legend()
plt.xlabel('EE80 segmented (arcsec)')
plt.ylabel('Quad. diff. (EE80 segmented - EE80 monolithic) / Wavelength')
plt.xlim([.1, .9])
plt.ylim([.0, 1e-3])
plt.savefig('results/ee80/ee80_quaddiff_wavscaled.png', bbox_inches='tight')
plt.close()
