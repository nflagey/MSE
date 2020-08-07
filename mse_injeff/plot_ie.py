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




def plotallvsiq():
    # - plot evolution as function of IQ for all configurations
    plt.figure()
    for zenith in zeniths:
        for field in fields:
            for specdiam in spectro_diams:
                for w in range(nwave):
                    ie_iqs = [ie_simu["iq" + str(iq)]['ZD' + zenith][field][specdiam][w] for iq in iqs]
                    plt.plot(iqs, ie_iqs / ie_iqs[0], color=wave_color[w], alpha=0.01)
    plt.title('IE as function of IQ')
    plt.xlabel('IQ (arcsec)')
    plt.ylabel('IE (relative)')
    plt.xlim([0.3, 2])
    plt.ylim([0, 1.25])
    plt.savefig('results/injeff_curves/simu_vs_IQ_all.png')
    plt.close()


def ploteachvsiq():
    # - plot evolution as function of IQ in one configuration, all wavelength
    for zenith in zeniths:
        for field in fields:
            for specdiam in spectro_diams:
                plt.figure()
                for w in range(nwave):
                    ie_iqs = [ie_simu["iq" + str(iq)]['ZD' + zenith][field][specdiam][w] for iq in iqs]
                    plt.plot(iqs, ie_iqs / ie_iqs[0], label=wave[w], color=wave_color[w])
                plt.title('IE as function of IQ')
                plt.xlabel('IQ (arcsec)')
                plt.ylabel('IE (relative)')
                plt.legend()
                plt.xlim([0.3, 2])
                plt.ylim([0, 1.25])
                plt.savefig('results/injeff_curves/simu_vs_IQ_' + field + specdiam + '_ZD' + zenith + '.png')
                plt.close()

#ieforbudget('optim')
ieforbudget('simu')
