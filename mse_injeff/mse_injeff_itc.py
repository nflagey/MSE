#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
import pickle

# Wavelengths
wave = np.asarray([360., 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784])  # in nm
nwave = len(wave)
# Weights for optimization
wave_weight = {'LR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 4, 2, 4, 2, 4, 1]),
               'MR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 3, 2, 1, 0, 0, 0]),
               'HR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 3, 2, 1, 0, 0, 0])}

# Define parameters for IQs
fields = ['X_+0.530', 'X_+0.650', 'X_+0.750', 'Y_+0.000',
          'Y_+0.375', 'Y_+0.530', 'Y_+0.650', 'Y_+0.750',
          'Y_-0.375', 'Y_-0.530', 'Y_-0.650', 'Y_-0.750']
zeniths = ['00', '30', '50', '60']
iqs = [0.3, 0.45, 0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.4, 1.55, 1.7, 1.85, 2.0]
spectro_diams = ['LR1.00', 'LR0.90', 'MR1.00', 'HR0.80', 'HR0.75', 'HR0.70']

# Positions in Zemax
yfields = ['Y_-0.750', 'Y_-0.650', 'Y_-0.530', 'Y_-0.375', 'Y_+0.000', 'Y_+0.375', 'Y_+0.530', 'Y_+0.650',
           'Y_+0.750']
xfields = ['X_+0.750', 'X_+0.650', 'X_+0.530', 'Y_+0.000', 'X_+0.530', 'X_+0.650', 'X_+0.750']
ypos = [-.75, -.65, -.53, -.375, 0, 0.375, 0.53, 0.65, 0.75]
xpos = [-.75, -.65, -.53, 0, 0.53, 0.65, 0.75]

# Open results
with open('results/no_segments_injeff_optim_curve_all.pkl', 'rb') as f:
    ie_optim = pickle.load(f)

with open('results/no_segments_injeff_curve_all.pkl', 'rb') as f:
    ie_simu = pickle.load(f)


def ieforitc():
    # Chose mean IE curves from simulations
    ie = ie_simu[0]  # [1] is for the standard deviation

    # Interpolate each IQ at 0.684º and find worst curve - these will be the curves used in the ITC
    # define fields
    ieposx, ieposy, ienegy = {}, {}, {}
    ieposx_wa, ieposy_wa, ienegy_wa = {}, {}, {}
    ieforitc = {}
    # loop over IQ
    for iq in iqs:
        striq = 'iq' + str(iq)
        ieposx[striq], ieposy[striq], ienegy[striq] = {}, {}, {}
        ieposx_wa[striq], ieposy_wa[striq], ienegy_wa[striq] = {}, {}, {}
        ieforitc[striq] = {}
        # loop over ZD
        for zenith in zeniths:
            strzd = 'ZD' + zenith
            ieposx[striq][strzd], ieposy[striq][strzd], ienegy[striq][strzd] = {}, {}, {}
            ieposx_wa[striq][strzd], ieposy_wa[striq][strzd], ienegy_wa[striq][strzd] = {}, {}, {}
            ieforitc[striq][strzd] = {}
            # loop over spectro
            for spectro in spectro_diams:
                ieposx[striq][strzd][spectro] = np.zeros_like(wave * 1.)
                ieposy[striq][strzd][spectro] = np.zeros_like(wave * 1.)
                ienegy[striq][strzd][spectro] = np.zeros_like(wave * 1.)
                # loop over wavelength
                for w in range(nwave):
                    # get IE(w) along X- and Y-axis
                    iealongx = np.asarray([ie[striq][strzd][xfield][spectro][w] for xfield in xfields])
                    iealongy = np.asarray([ie[striq][strzd][yfield][spectro][w] for yfield in yfields])
                    # interpol at 0.684º
                    ieposx[striq][strzd][spectro][w] = np.asarray([np.interp(0.684, xpos, iealongx)])
                    ieposy[striq][strzd][spectro][w] = np.asarray([np.interp(0.684, ypos, iealongy)])
                    ienegy[striq][strzd][spectro][w] = np.asarray([np.interp(-0.684, ypos, iealongy)])

                # weighted average
                ieposx_wa[striq][strzd][spectro] = np.sum(ieposx[striq][strzd][spectro] * wave_weight[spectro[0:2]])\
                                                   / np.sum(wave_weight[spectro[0:2]])
                ieposy_wa[striq][strzd][spectro] = np.sum(ieposy[striq][strzd][spectro] * wave_weight[spectro[0:2]])\
                                                   / np.sum(wave_weight[spectro[0:2]])
                ienegy_wa[striq][strzd][spectro] = np.sum(ienegy[striq][strzd][spectro] * wave_weight[spectro[0:2]])\
                                                   / np.sum(wave_weight[spectro[0:2]])
                # pick the worst one for IE budget
                argmin = np.argmin([ieposx_wa[striq][strzd][spectro], ieposy_wa[striq][strzd][spectro],
                                    ienegy_wa[striq][strzd][spectro]])
                ieforitc[striq][strzd][spectro] = [ieposx[striq][strzd][spectro], ieposy[striq][strzd][spectro],
                                                   ienegy[striq][strzd][spectro]][argmin]

    # Save as Numpy file to avoid pickle protocol issue in the ITC
    np.save('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/no_segments_injeff_curve_for_itc_all.npy', ieforitc)


def ieforbudget(type='simu'):
    # Chose IE simu or IE optim
    if type == 'simu':
        ie = ie_simu[0]  # [1] is for the standard deviation
    else:
        if type == 'sdev':
            ie = ie_simu[1]  # [1] is for the standard deviation
        else:
            ie = ie_optim

    # Get IE at each position at 0.56" IQ budget
    ie056 = {}
    for spectro in spectro_diams:
        ie056[spectro] = {}
        print(spectro)
        for field in fields:
            ie045 = [float(ie['iq0.45']['ZD30'][field][spectro][w]) for w in range(nwave)]
            ie060 = [float(ie['iq0.6']['ZD30'][field][spectro][w]) for w in range(nwave)]
            ie056[spectro][field] = [np.interp(0.56, [0.45, 0.6], [ie045[w], ie060[w]]) for w in range(nwave)]
            print(field, ie056[spectro][field])

    # Interpolate at 0.684º
    # define fields
    ie056posx, ie056posy, ie056negy = {}, {}, {}
    ie056posx_wa, ie056posy_wa, ie056negy_wa = {}, {}, {}
    ie056budget = {}
    # loop over spectro
    for spectro in spectro_diams:
        # interpol at 0.684º
        ie056posx[spectro] = np.asarray([np.interp(0.684, xpos, np.asarray([ie056[spectro][xfield][w]
                                                                            for xfield in xfields]))
                                         for w in range(nwave)])
        ie056posy[spectro] = np.asarray([np.interp(0.684, ypos, np.asarray([ie056[spectro][yfield][w]
                                                                            for yfield in yfields]))
                                         for w in range(nwave)])
        ie056negy[spectro] = np.asarray([np.interp(-0.684, ypos, np.asarray([ie056[spectro][yfield][w]
                                                                             for yfield in yfields]))
                                         for w in range(nwave)])
        # weighted average
        ie056posx_wa[spectro] = np.sum(ie056posx[spectro] * wave_weight[spectro[0:2]]) / np.sum(wave_weight[spectro[0:2]])
        ie056posy_wa[spectro] = np.sum(ie056posy[spectro] * wave_weight[spectro[0:2]]) / np.sum(wave_weight[spectro[0:2]])
        ie056negy_wa[spectro] = np.sum(ie056negy[spectro] * wave_weight[spectro[0:2]]) / np.sum(wave_weight[spectro[0:2]])
        # pick the worst one for IE budget
        ie056budget[spectro] = [ie056posx[spectro], ie056posy[spectro], ie056negy[spectro]][np.argmax([ie056posx_wa[spectro], ie056posy_wa[spectro], ie056negy_wa[spectro]])]

    # print at zemax wavelengths
        print(spectro, ie056budget[spectro])

    # interpol at other budget wavelengths
    budget_waves = np.asarray([360, 370, 400, 482, 626, 767, 900, 910, 950, 962, 1235, 1300, 1500, 1662, 1800])
    ie056budgetwaves = {}
    # loop over wavelength
    for spectro in spectro_diams:
        ie056budgetwaves[spectro] = np.interp(budget_waves, wave, ie056budget[spectro])
        print(spectro, ie056budgetwaves[spectro])

    print('hello')


# Execute function
ieforitc()
ieforbudget('simu')