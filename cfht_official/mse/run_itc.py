#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
from mse_itc_dev import MseSpectrum
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

path = '/mse_itc_dev/plots/'

# Create spectrum
spectro = 'LR'
# -- LR
if spectro == 'LR':
    s0 = MseSpectrum(tgtmag=24, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                     coating='ZeCoat', fibdiam=1.0, spectro='LR', spatbin=2, specbin=1, meth='getSNR', snr=1, etime=3600,
                     badtgtmag=19)
# -- MR
elif spectro == 'MR':
    s0 = MseSpectrum(tgtmag=23.5, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.56,
                     coating='ZeCoat', fibdiam=1.0, spectro='MR', spatbin=2, specbin=1, meth='getSNR', snr=1, etime=3600,
                     badtgtmag=19)
# -- HR
else:
    s0 = MseSpectrum(tgtmag=20, band='g', template='flat', redshift=0, airmass=1.2, skymag=19.5, seeing=0.56,
                     coating='ZeCoat', fibdiam=0.8, spectro='HR', spatbin=2, specbin=1, meth='getSNR', snr=1, etime=3600,
                     badtgtmag=16)

# Check intermediate steps
# Apply extinction
s1 = s0.apply_atmos_ext()
# Apply frontend throughput
s2 = s1.apply_throughput_front()
#plt.semilogy(s2.wgrid, s2.skyflux)
#plt.plot(s2.wgrid, s2.tgtflux)
#plt.title("At the focal surface")
#plt.show()

# Apply injection efficiency
s3 = s2.apply_injeff()
# Apply backend throughput
s4 = s3.apply_throughput_back()

# Compute SNR at once
s5, script, div = s0.compute_snr(doplot='offline')


# Verify numbers are consistent with budgets
np.set_printoptions(precision=1)

lam1 = np.array([360, 370, 400, 482, 626, 767, 900, 910, 950, 962, 1235, 1300, 1500, 1662, 1800]) * 10.

plt.figure()
for i in range(1 + int(np.max(s5.armgrid))):
    arm = s5.armgrid == i
    # Total is quadrature sum of all sources of noise
    plt.semilogy(s5.wgrid[arm], signal.medfilt(np.sqrt(signal.medfilt(s5.skynoise[arm], 101)**2
                                                       + signal.medfilt(s5.tgtnoise[arm], 101)**2
                                                       + s5.xtalk[arm] + s5.ghost[arm] + s5.darknoise[arm]**2
                                                       + s5.teldiffuse[arm] + s5.instdiffuse[arm]
                                                       + s5.thermalnoise[arm]**2 + s5.readout[arm]**2), 101),
                 label='Total', c='#000000')

    plt.semilogy(s5.wgrid[arm], signal.medfilt(s5.skynoise[arm], 101), label='Sky', c='#0000FF')

    plt.semilogy(s5.wgrid[arm], signal.medfilt(s5.tgtnoise[arm], 101), label='Target', c='#FF0000')

    plt.semilogy(s5.wgrid[arm], signal.medfilt(np.sqrt(s5.xtalk[arm]), 101), label='X-talk', linestyle='dashed', c='#FF6600')

    plt.semilogy(s5.wgrid[arm], signal.medfilt(np.sqrt(s5.ghost[arm]), 101), label='Ghosts', linestyle='dotted', c='#0066FF', )

    plt.semilogy(s5.wgrid[arm], s5.darknoise[arm], label='Dark', c='#880088')

    plt.semilogy(s5.wgrid[arm], np.sqrt(s5.teldiffuse[arm]), label='Diff. tel.', c='#444444')

    plt.semilogy(s5.wgrid[arm], np.sqrt(s5.instdiffuse[arm]), label='Diff. instr.', c='#888888')

    plt.semilogy(s5.wgrid[arm], s5.thermalnoise[arm], label='Thermal', c='#BB8800')

    plt.semilogy(s5.wgrid[arm], s5.readout[arm], label='Readout', c='#00FF00')

    if i == 0:
        plt.legend(loc="lower right", bbox_to_anchor=(1.25, 0))

plt.xlabel('Wavelength (A)')
plt.ylabel('Counts (e/res.elem.)')
plt.title("Noise")
plt.savefig(path + 'noise_'+s0.spectro+'.png', bbox_inches="tight")

plt.figure()
for i in range(1 + int(np.max(s4.armgrid))):
    arm = s4.armgrid == i
    plt.plot(s5.wgrid[arm], s5.snr[arm], 'k', alpha=0.25)
    plt.plot(s5.wgrid[arm], signal.medfilt(s5.snr[arm], 101), 'k')
if spectro == 'LR':
    plt.plot([3600, 3700, 4000, 4000, 18000], [1,1,1,2,2], 'c:')
elif spectro == 'MR':
    plt.plot([3600, 3700, 4000, 4000, 9500], [1, 1, 1, 2, 2], 'c:')
else:
    plt.plot([3600, 3700, 4000, 4000, 9000], [5, 5, 5, 10, 10], 'c:')
plt.title("SNR")
plt.savefig(path + 'snr'+s0.spectro+'.png', bbox_inches="tight")


print('Done')