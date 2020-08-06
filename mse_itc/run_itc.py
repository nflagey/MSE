#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
from mse_itc import MseSpectrum
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from astropy.table import Table

path = 'output/'

# IE wavelength grid
lam0 = np.array([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784]) * 10.
# Sensitivity wavelength grid
lam1 = np.array([360, 370, 400, 482, 626, 767, 900, 910, 950, 962, 1235, 1300, 1500, 1662, 1800]) * 10.

# Create spectrum
spectro = 'LR'
# -- LR
if spectro == 'LR':
    spec = MseSpectrum(tgtmag=24, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                       coating='ZeCoat', fibdiam=1.0, spectro='LR', spatbin=2, specbin=1, meth='getSNR', snr=1,
                       etime=3600, badtgtmag=19)
# -- MR
elif spectro == 'MR':
    spec = MseSpectrum(tgtmag=23.5, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                       coating='ZeCoat', fibdiam=1.0, spectro='MR', spatbin=2, specbin=1, meth='getSNR', snr=1,
                       etime=3600, badtgtmag=19)
# -- HR
else:
    spec = MseSpectrum(tgtmag=20, band='g', template='flat', redshift=0, airmass=1.2, skymag=19.5, seeing=0.56,
                       coating='ZeCoat', fibdiam=0.8, spectro='HR', spatbin=2, specbin=1, meth='getSNR', snr=1,
                       etime=3600, badtgtmag=16)

# -- Compute SNR at once (compute_snr calls all the previous steps)
script, div = spec.compute_snr(doplot='offline')  # use doplot='no' to avoid creating the html page with plots

# -- Some plots and tables
# Table with sky counts at focal surface
t = Table((spec.wgrid / 10, spec.tgtflux['at_focus']), names=('wave', 'targetflux'))
t.write(path + 'target_flux.dat', format='ascii', overwrite=True)

t = Table((spec.wgrid, spec.skyflux['at_focus'] / spec.egrid, spec.tgtflux['at_focus'] / spec.egrid),
          names=('wave', 'skycount', 'targetcount'))
t.write(path + 'counts_at_telescope_focal_surface.dat', format='ascii', overwrite=True)

# Plot counts at detector
# plt.semilogy(spec.wgrid, spec.skydetec/spec.etime, label='Sky')
# plt.plot(spec.wgrid, spec.tgtdetec/spec.etime, label='Target')
# plt.title("At the detector")
# plt.ylabel("Photons/s per resolution element")
# plt.legend()
# plt.show()


# -- Verify numbers are consistent with budgets and main plot
np.set_printoptions(precision=1)

plt.figure()
for i in range(1 + int(np.max(spec.armgrid))):
    arm = spec.armgrid == i
    # Total is quadrature sum of all sources of noise
    plt.semilogy(spec.wgrid[arm], signal.medfilt(np.sqrt(signal.medfilt(spec.skynoise[arm], 101)**2
                                                         + signal.medfilt(spec.tgtnoise[arm], 101)**2
                                                         + spec.xtalk[arm] + spec.ghost[arm] + spec.darknoise[arm]**2
                                                         + spec.teldiffuse[arm] + spec.instdiffuse[arm]
                                                         + spec.thermalnoise[arm]**2 + spec.readout[arm]**2), 101),
                 label='Total', c='#000000')

    plt.semilogy(spec.wgrid[arm], signal.medfilt(spec.skynoise[arm], 101), label='Sky', c='#0000FF')
    plt.semilogy(spec.wgrid[arm], signal.medfilt(spec.tgtnoise[arm], 101), label='Target', c='#FF0000')
    plt.semilogy(spec.wgrid[arm], signal.medfilt(np.sqrt(spec.xtalk[arm]), 101), label='X-talk', linestyle='dashed', c='#FF6600')
    plt.semilogy(spec.wgrid[arm], signal.medfilt(np.sqrt(spec.ghost[arm]), 101), label='Ghosts', linestyle='dotted', c='#0066FF', )
    plt.semilogy(spec.wgrid[arm], spec.darknoise[arm], label='Dark', c='#880088')
    plt.semilogy(spec.wgrid[arm], np.sqrt(spec.teldiffuse[arm]), label='Diff. tel.', c='#444444')
    plt.semilogy(spec.wgrid[arm], np.sqrt(spec.instdiffuse[arm]), label='Diff. instr.', c='#888888')
    plt.semilogy(spec.wgrid[arm], spec.thermalnoise[arm], label='Thermal', c='#BB8800')
    plt.semilogy(spec.wgrid[arm], spec.readout[arm], label='Readout', c='#00FF00')

    # Print values for budget
    print("Arm", i)
    print("Dark", np.interp(lam1, spec.wgrid[arm], spec.dark[arm]))
    print("Readout", np.interp(lam1, spec.wgrid[arm], spec.readout[arm]))

    if i == 0:
        plt.legend(loc="lower right", bbox_to_anchor=(1.25, 0))

plt.xlabel('Wavelength (A)')
plt.ylabel('Counts (e/res.elem.)')
plt.title("Noise")
plt.savefig(path + 'noise_'+spec.spectro+'.png', bbox_inches="tight")

plt.figure()
for i in range(1 + int(np.max(spec.armgrid))):
    arm = spec.armgrid == i
    plt.plot(spec.wgrid[arm], spec.snr[arm], 'k', alpha=0.25)
    plt.plot(spec.wgrid[arm], signal.medfilt(spec.snr[arm], 101), 'k')
if spectro == 'LR':
    plt.plot([3600, 3700, 4000, 4000, 18000], [1, 1, 1, 2, 2], 'c:')
elif spectro == 'MR':
    plt.plot([3600, 3700, 4000, 4000, 9500], [1, 1, 1, 2, 2], 'c:')
else:
    plt.plot([3600, 3700, 4000, 4000, 9000], [5, 5, 5, 10, 10], 'c:')
plt.title("SNR")
plt.savefig(path + 'snr_'+spec.spectro+'.png', bbox_inches="tight")


print('Done')