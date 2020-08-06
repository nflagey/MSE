#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
from mse_itc import MseSpectrum
import numpy as np
from scipy import signal
from bokeh.plotting import figure, show

# Spatial binning
spatbin = 2

# Create spectrum according to LR sensitivity requirement for CoDR design
spec_codr = MseSpectrum(tgtmag=24, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                        coating='ZeCoat', fibdiam=1.0, spectro='LR', spatbin=spatbin, specbin=1, meth='getSNR', snr=1,
                        etime=3600, badtgtmag=19, lmr='codr')
# Compute SNR at once
script_codr, div_codr = spec_codr.compute_snr(doplot='online')
# Plot
fig = figure(title="Spectra", y_axis_label="Flux (erg/s/cm2/A)", x_axis_label="Wavelength (A)",
             x_range=(3500, 18000), y_range=(0, 5))
for i in range(1 + int(np.max(spec_codr.armgrid))):
    arm = spec_codr.armgrid == i
    fig.line(spec_codr.wgrid[arm], signal.medfilt(spec_codr.snr[arm], 101), line_color='blue', legend_label='CoDR - LR')

# Create spectrum according to LR sensitivity requirement for "MOONS" design
spec_moons = MseSpectrum(tgtmag=24., band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                         coating='ZeCoat', fibdiam=1.0, spectro='LR', spatbin=spatbin, specbin=1, meth='getSNR', snr=1,
                         etime=3600, badtgtmag=19, lmr='moons')
# Compute SNR at once
script_moons, div_moons = spec_moons.compute_snr(doplot='online')
# Plot
for i in range(1 + int(np.max(spec_moons.armgrid))):
    arm = spec_moons.armgrid == i
    fig.line(spec_moons.wgrid[arm], signal.medfilt(spec_moons.snr[arm], 101), line_color='#AA0000', legend_label='MOONS-like')

# Bin result
for i in range(1 + int(np.max(spec_moons.armgrid))):
    arm = spec_moons.armgrid == i
    spec_bin_wav = np.array([0.5*(spec_moons.wgrid[arm][j] + spec_moons.wgrid[arm][j+1])
                             for j in range(0, len(spec_moons.wgrid[arm])-1, 2)])
    spec_bin_snr = np.array([(spec_moons.snr[arm][j] + spec_moons.snr[arm][j+1])/np.sqrt(2)
                             for j in range(0, len(spec_moons.wgrid[arm])-1, 2)])
# Plot
    fig.line(spec_bin_wav, signal.medfilt(spec_bin_snr, 101), line_color='#FF0000', legend_label='MOONS-like binned')


fig.axis.major_label_text_font_size = "12pt"
fig.axis.axis_label_text_font_size = "14pt"


show(fig)

print('Done!')