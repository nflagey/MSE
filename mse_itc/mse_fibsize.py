#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
# from mse_itc.mse_itc import *
from mse_itc import *
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
import pickle
from astropy.table import Table
import os

# Need to set the current path in case this code is called from somewhere else
here_path = os.path.dirname(os.path.realpath(__file__))

# === PARAMETERS ===
# Magnitudes
tgt_mag = {'LR': 24, 'MR': 23.5, 'HR': 20}
sky_mag = {'LR': 20.7, 'MR': 20.7, 'HR': 19.5}
# Airmass
am_vals = {'ZD00': 1.0, 'ZD30': 1.2, 'ZD50': 1.5}  # TODO those are approximate and should be revised at some point


# === Compute and save individual SNR curves ===
def comp_snr(iq_vals=None, zd_vals=None, fibers=None, src_type='point'):
    # Initiate SNR dictionary
    res = {}
    # Loop over IQs
    for iq_val in iq_vals:
        res[iq_val] = {}
        # Loop over ZD
        for zd_val in zd_vals:
            res[iq_val][zd_val] = {}
            plt.figure()
            # Loop over fiber sizes
            for fiber in fibers:
                # - create spectrum
                spectro = fiber[0:2]  # LR, MR, HR ... should be the same for all
                fib_size = fiber[2:]  # 1.00, 0.75, ... as a string
                spec = MseSpectrum(sessionID=-1, tgtmag=tgt_mag[spectro], band='g', template='flat',
                                   src_type=src_type, redshift=0, airmass=am_vals[zd_val], skymag=sky_mag[spectro],
                                   seeing=float(iq_val[2:]), coating='ZeCoat', fibdiam=float(fib_size),
                                   spectro=spectro, specbin=1, spatbin=1, meth='getSNR', snr=1, etime=3600)
                # - compute SNR
                script, div = spec.compute_snr(doplot='online')
                # - store into dictionary
                res[iq_val][zd_val][fiber] = spec
                # - plot SNR ratio of current fiber to first fiber size
                if fiber != fibers[0]:
                    plt.plot(spec.wgrid, res[iq_val][zd_val][fiber].snr / res[iq_val][zd_val][fibers[0]].snr,
                             label=iq_val + ' ' + zd_val + ' ' + fib_size, alpha=0.5)
            plt.ylabel('SNR ratio')
            plt.ylim([0.9, 1.1])
            plt.legend()
            plt.savefig(here_path + '/data_fibsize/SNR_' + iq_val + '_' + zd_val + '.png')

    # Save individual SNR curves
    with open(here_path + '/data_fibsize/fibsize_snrcomp_' + spectro + '_' + src_type + '.pkl', 'wb') as output:
        pickle.dump(res, output, pickle.HIGHEST_PROTOCOL)


# === Plot ratio for all IQ ===
# ... at each ZD independently, one fiber size at a time
def plot_snr_ratio_all(iq_vals=None, zd_vals=None, fibers=None, src_type='point'):
    # Restore results
    spectro = fibers[0][0:2]
    with open(here_path + '/data_fibsize/fibsize_snrcomp_' + spectro + '_' + src_type + '.pkl', 'rb') as input:
        res = pickle.load(input)
    # Color variable
    color = np.asarray(range(len(iq_vals))) / len(iq_vals)
    color += 0.5 * color[1]
    # Loop over fiber size (all but the first, since we do ratios)
    for fiber in fibers[1:]:
        fig = plt.figure(figsize=(13.66, 7.68), tight_layout=True)
        count = 1
        # Loop over zenith distance
        for zd_val in zd_vals:
            ax1 = fig.add_subplot(2, 2, count)
            # Loop over IQs
            for i, iq_val in enumerate(iq_vals):
                # Get results for both fiber size
                snr = res[iq_val][zd_val][fiber].snr
                snr0 = res[iq_val][zd_val][fibers[0]].snr
                wav = res[iq_val][zd_val][fibers[0]].wgrid
                armgrid = res[iq_val][zd_val][fibers[0]].armgrid
                # Plot
                for j in range(1 + int(np.max(armgrid))):
                    arm = (armgrid == j)
                    ax1.plot(wav[arm], (snr / snr0)[arm], c=(color[i], 0.5, 0.5), alpha=0.1)
                    if j == 1 and count == 3:
                        ax1.plot(wav[arm], signal.medfilt((snr / snr0)[arm], 101), c=(color[i], 0.5, 0.5), label=iq_val)
                    else:
                        ax1.plot(wav[arm], signal.medfilt((snr / snr0)[arm], 101), c=(color[i], 0.5, 0.5))
            ax1.set_title(zd_val)
            ax1.set_xlabel('Wavelength (A)')
            ax1.set_ylabel('SNR ratio ('+fiber+'" to '+fibers[0]+'")')
            ax1.set_ylim([0.95, 1.25])
            if count == 3:
                ax1.legend(fontsize=7, bbox_to_anchor=(1.4, 1.0))
            count += 1
        fig.savefig(here_path + '/data_fibsize/fibsize_snrcomp' + spectro + '_' + fiber + 'to' + fibers[0] + '.png')


# === Fold in the IQ distribution ===
def plot_snr_ratio_iq(iq_vals=None, zd_vals=None, fibers=None, src_type='point'):
    # Restore results
    spectro = fibers[0][0:2]
    with open(here_path + '/data_fibsize/fibsize_snrcomp_' + spectro + '_' + src_type + '.pkl', 'rb') as input:
        res = pickle.load(input)
    plt.figure()
    # - get IQ and AM values
    mkam = Table.read(here_path + '/data_fibsize/mkam_dimm_all.csv', format='csv')
    mkam_iq_in = mkam['massdimm_seeing']
    mkam_am_in = mkam['massdimm_airmass']
    # -- remove unwanted values
    mkam_iq_out = [float(mkam_iq_in[i]) for i in range(len(mkam_iq_in)) if mkam_iq_in[i] != '#NAME?']
    mkam_am_out = [float(mkam_am_in[i]) for i in range(len(mkam_iq_in)) if mkam_iq_in[i] != '#NAME?']
    # -- convert into numpy.arrays
    iq = np.asarray(mkam_iq_out)
    am = np.asarray(mkam_am_out)
    # - get IQ_500
    # -- account for airmass to go back to zenith
    iq = iq / am**(3/5)
    # -- remove GL
    iq = (iq**(5/3) - 0.289**(5/3))**(3/5)
    # -- compute r0 at 500nm, with IQ in rad
    r0 = 0.98 * 500e-9 / (iq / 180 * np.pi / 3600)
    # -- convert MKAM-DIMM values into IQ_500 values with L0 = 30m
    iq00 = iq * np.sqrt(1 - 2.183 * (r0 / 30)**0.356)
    # - account for airmass to go to ZD30 and ZD50
    am30 = 1 / np.cos(30 / 180 * np.pi)
    am50 = 1 / np.cos(50 / 180 * np.pi)
    iq30 = iq * am30**(3/5)
    iq50 = iq * am50**(3/5)
    # - print median values
    print(np.nanmedian(iq00), np.nanmedian(iq30), np.nanmedian(iq50))
    # - get MSE IQ at focal surface (that used for IE calculations) by adding GL uplift and Thermal
    iq00 = (iq00 ** (5./3) + 0.2 ** (5./3) + 0.1 ** (5./3)) ** (3./5)
    iq30 = (iq30 ** (5./3) + 0.2 ** (5./3) + 0.1 ** (5./3)) ** (3./5)
    iq50 = (iq50 ** (5./3) + 0.2 ** (5./3) + 0.1 ** (5./3)) ** (3./5)
    # - make distributions at the three ZD
    h00 = np.histogram(iq00, bins=21, range=(-0.025, 1.025))  # 88% of the IQ values are below 1.025"
    h30 = np.histogram(iq30, bins=21, range=(-0.025, 1.025))
    h50 = np.histogram(iq50, bins=21, range=(-0.025, 1.025))
    norm00 = np.sum(h00[0])
    norm30 = np.sum(h30[0])
    norm50 = np.sum(h50[0])

    # Plot those distributions
    fig = plt.figure(figsize=(13.66, 7.68), tight_layout=True)
    plt.hist(iq00, bins=21, range=(-0.025, 1.025), color='#880000', label='ZD00', fill=False, histtype='step')
    plt.hist(iq00, bins=41, range=(-0.025, 2.025), color='#FF0000', alpha=0.5, fill=False, histtype='step')
    plt.hist(iq30, bins=21, range=(-0.025, 1.025), color='#008800', label='ZD30', fill=False, histtype='step')
    plt.hist(iq30, bins=41, range=(-0.025, 2.025), color='#00FF00', alpha=0.5, fill=False, histtype='step')
    plt.hist(iq50, bins=21, range=(-0.025, 1.025), color='#000088', label='ZD50', fill=False, histtype='step')
    plt.hist(iq50, bins=41, range=(-0.025, 2.025), color='#0000FF', alpha=0.5, fill=False, histtype='step')
    plt.title('Distribution of IQ at MSE focal surface')
    plt.legend()
    fig.savefig(here_path + '/data_fibsize/fibsize_snrcomp_iqdist.png')

    # Loop to build SNR over IQ distribution
    wav = res[iq_vals[0]][zd_vals[0]][fibers[0]].wgrid
    meansnr = {'ZD00': {}, 'ZD30': {}, 'ZD50': {}}
    for fiber in fibers:
        for zd_val in zd_vals[0:3]:
            meansnr[zd_val][fiber] = np.zeros_like(wav)
        for i in range(len(iq_vals)):
            iq_val = iq_vals[i]
            meansnr['ZD00'][fiber] += res[iq_val]['ZD00'][fiber].snr * h00[0][i] / norm00
            meansnr['ZD30'][fiber] += res[iq_val]['ZD30'][fiber].snr * h30[0][i] / norm30
            meansnr['ZD50'][fiber] += res[iq_val]['ZD50'][fiber].snr * h50[0][i] / norm50

        table = Table([wav, meansnr['ZD00'][fiber]], names=('wav', 'snr'))
        table.write(here_path + '/data_fibsize/fibsize_snrcomp_iqdist' + spectro + '_' + fiber + '_ZD00.txt',
                    format='ascii')
        table = Table([wav, meansnr['ZD30'][fiber]], names=('wav', 'snr'))
        table.write(here_path + '/data_fibsize/fibsize_snrcomp_iqdist' + spectro + '_' + fiber + '_ZD30.txt',
                    format='ascii')
        table = Table([wav, meansnr['ZD50'][fiber]], names=('wav', 'snr'))
        table.write(here_path + '/data_fibsize/fibsize_snrcomp_iqdist' + spectro + '_' + fiber + '_ZD50.txt',
                    format='ascii')

    # Plot of ratio, at three different ZD
    armgrid = res[iq_vals[0]][zd_vals[0]][fibers[0]].armgrid
    # Color variable
    color = [0.1, 0.5, 0.9]
    # loop over fiber size
    for fiber in fibers[1:]:
        fig = plt.figure(figsize=(13.66, 7.68), tight_layout=True)
        count = 0
        for zd_val in zd_vals[0:3]:
            for j in range(1 + int(np.max(armgrid))):
                arm = (armgrid == j)
                if j == 0:
                    plt.plot(wav[arm], (meansnr[zd_val][fiber])[arm] / (meansnr[zd_val][fibers[0]])[arm],
                             c=(color[count], 0.5, 0.5), label=zd_val, alpha=0.5)
                    plt.plot(wav[arm],
                             signal.medfilt((meansnr[zd_val][fiber])[arm] / (meansnr[zd_val][fibers[0]])[arm], 101),
                             c=(color[count], 0.5, 0.5))
                else:
                    plt.plot(wav[arm], (meansnr[zd_val][fiber])[arm] / (meansnr[zd_val][fibers[0]])[arm],
                             c=(color[count], 0.5, 0.5), alpha=0.5)
                    plt.plot(wav[arm],
                             signal.medfilt((meansnr[zd_val][fiber])[arm] / (meansnr[zd_val][fibers[0]])[arm], 101),
                             c=(color[count], 0.5, 0.5))
                plt.title(fiber+' / '+fibers[0])
                plt.legend()
                plt.xlabel('Wavelength (A)')
                plt.ylabel('SNR ratio (' + fiber + '" to ' + fibers[0] + '")')
                plt.ylim([0.95, 1.25])
            count += 1

        fig.savefig(here_path + '/data_fibsize/fibsize_snrcomp_iqdist' + spectro + '_' + fiber + 'to' + fibers[0] + '.png')
