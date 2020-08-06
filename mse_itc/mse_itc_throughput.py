#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
from astropy.table import Table

# Path for Throughput files
path = 'THROUGHPUT/'

# Wavelengths in Angstroms (3500 - 18500)
lam = np.arange(1500, dtype=float) * 10 + 3500.
# IE wavelength grid
lam0 = np.array([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784]) * 10.
# Sensitivity wavelength grid
lam1 = np.array([360, 370, 400, 482, 626, 767, 900, 910, 950, 962, 1235, 1300, 1500, 1662, 1800]) * 10.


def mse_itc_throughput():
    """Writes the throughput files for the ITC
    """

    # = = = = = = = = =
    # Structure
    # = = = = = = = = =

    struc_thr = np.ones_like(lam, dtype=float) * 0.96
    # save file
    struc_data = Table([lam, struc_thr], names=['lamA', 'thr'])
    ascii.write(struc_data, path + 'mse_itc_throughput_struc.dat', overwrite=True)
    # print values at Lam1 for budget
    print('Structure', np.interp(lam1, lam, struc_thr))

    # = = = = = = = = =
    # M1 (ZeCoat and Gemini curves, then remove 1.0% to account for time until recoating and 0.6% for edge bevel)
    # = = = = = = = = =

    # ZeCoat
    m1 = Table.read(path + 'zecoat_zc562.dat', format='ascii')
    m1_thr = np.interp(lam, m1['col1'].data * 10, m1['col2'].data) / 100
    m1_thr *= .99 * .994
    # save file
    m1_data = Table([lam, m1_thr], names=['lamA', 'thr'])
    ascii.write(m1_data, path + 'mse_itc_throughput_m1_zecoat.dat', overwrite=True)
    # print values at Lam0 and Lam1 for budget
    print('M1 ZeCoat', np.interp(lam1, lam, m1_thr))

    # Gemini
    m1_gem = Table.read(path + 'Gemini_reflectivity.dat', format='ascii')
    m1_gem_thr = np.interp(lam, m1_gem['col1'].data * 10,
                       m1_gem['col2'].data + m1_gem['col3'].data + m1_gem['col4'].data + m1_gem['col5'].data) / 400
    m1_gem_thr *= .99 * .994
    # save file
    m1_gem_data = Table([lam, m1_gem_thr], names=['lamA', 'thr'])
    ascii.write(m1_gem_data, path + 'mse_itc_throughput_m1_gemini.dat', overwrite=True)
    # print values at Lam0 and Lam1 for budget
    print('M1 Gemini', np.interp(lam1, lam, m1_gem_thr))

    # plot and save
    plt.figure(figsize=(10, 5))
    plt.semilogx(lam/10., m1_thr, label='ZeCoat')
    plt.plot(lam/10., m1_gem_thr, label='Gemini')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Reflectivity')
    plt.xlim([360, 1800])
    plt.ylim([0.50, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_m1_zecoat.png')

    # = = = = = = = = =
    # Prime focus unit
    # = = = = = = = = =

    # AR coating, degrading it by 20% for contingency
    arsili = Table.read(path + 'solgel_arcoating_sili.csv', format='csv')
    sili_arcoat = np.interp(lam, arsili['wave'].data * 10, 1 - arsili['reflec'].data / 100 * 1.2)
    arpbm = Table.read(path + 'solgel_arcoating_pbm.csv', format='csv')
    pbm_arcoat = np.interp(lam, arpbm['wave'].data * 10, 1 - arpbm['reflec'].data / 100 * 1.2)

    # obscuration
    pfue_obsc = np.ones_like(lam) * 0.987

    # vignetting at 90% radius
    pfue_vign = np.ones_like(lam) * 0.91

    # absorption
    # PBM2Y from Ohara
    pfue_wav_pb = np.array([360., 365., 370., 380., 390., 400., 420., 440., 460., 480., 500.,
                            550., 600., 650., 700., 800., 900., 1000, 1200, 1400, 1600, 1800])
    pfue_trans_pb_25 = np.array([.951, .965, .978, .987, .991, .993, .995, .995, .996, .996, .997,
                                 .998, .998, .997, .998, .998, .998, .995, .995, .990, .985, .951])
    pfue_pb_thick = 100.  # mm
    pfue_abs_pb = (1. - pfue_trans_pb_25) / 25.  # per mm
    pfue_trans_pb = np.exp(-pfue_abs_pb * pfue_pb_thick)  # for total thickness
    pfue_trans_pb = np.interp(lam, pfue_wav_pb * 10., pfue_trans_pb)  # transmission
    #  Fused silica from Corning 7980
    pfue_wav_sili = np.array([360., 365., 370., 380., 390., 400., 420., 440., 460., 480.,
                              500., 550., 600., 650., 700., 800., 900., 1000, 1200, 1210,
                              1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310,
                              1320, 1330, 1340, 1350, 1360, 1370, 1380, 1390, 1400, 1410,
                              1420, 1430, 1440, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800])
    pfue_trans_sili_10 = np.array([.9995, .9995, .9997, .9997, .9998, .9999, .9999, .9999, .9998, .9997,
                                   .9996, .9994, .9994, .9994, .9992, .9991, .9990, .9989, .9982, .9980,
                                   .9977, .9965, .9930, .9930, .9950, .9963, .9967, .9969, .9968, .9965,
                                   .9956, .9937, .9908, .9881, .9846, .9575, .8824, .9040, .9450, .9685,
                                   .9802, .9886, .9908, .9934, .9972, .9975, .9974, .9973, .9970, .9965, .9958])
    pfue_sili_thick = 300.  # mm
    pfue_abs_sili = (1. - pfue_trans_sili_10) / 10.  # per mm
    pfue_trans_sili = np.exp(-pfue_abs_sili * pfue_sili_thick)  # for total thickness
    pfue_trans_sili = np.interp(lam, pfue_wav_sili * 10., pfue_trans_sili)  # transmission

    # Total for PFUE
    pfue_sili_nlens = 3
    pfue_pb_nlens = 2
    # -- for ITC: without vignetting because it is accounted for in the IE curves ---> not true anymore!!!
    #    though need to remove obscuration due to PFHS, to be consistent with Budget
    pfue_thr_itc = sili_arcoat ** (pfue_sili_nlens * 2) * pbm_arcoat ** (pfue_pb_nlens * 2) *\
                   pfue_trans_sili * pfue_trans_pb * pfue_vign
    # save file for ITC
    pfue_data = Table([lam, pfue_thr_itc], names=['lamA', 'thr'])
    ascii.write(pfue_data, path + 'mse_itc_throughput_pfue.dat', overwrite=True)
    pfhs_data = Table([lam, pfue_obsc], names=['lamA', 'thr'])
    ascii.write(pfhs_data, path + 'mse_itc_throughput_pfhs.dat', overwrite=True)
    # -- for budget: with vignetting
    pfue_thr_bud = sili_arcoat ** (pfue_sili_nlens * 2) * pbm_arcoat ** (pfue_pb_nlens * 2) *\
                   pfue_trans_sili * pfue_trans_pb * pfue_vign
    # print values at Lam1 for budget
    print('WFC/ADC reflec.', np.interp(lam1, lam, sili_arcoat ** (pfue_sili_nlens * 2) * pbm_arcoat ** (pfue_pb_nlens * 2)))
    print('WFC/ADC trans.', np.interp(lam1, lam, pfue_trans_sili * pfue_trans_pb))
    print('WFC/ADC vign.', np.interp(lam1, lam, pfue_vign))
    # obscuration is part of PFHS in the budget --- but vignetting is part of WFC/ADC in the budget
    print('PFHS obsc.', np.interp(lam1, lam, pfue_obsc))
    print('WFC/ADC total', np.interp(lam1, lam, pfue_thr_bud))

    # Read throughput file obtained from Will's spreadsheet
    pfue_ws = ascii.read(path + 'mse_itc_throughput_pfue_WS.dat', format= 'csv')
    # print values at Lam1 for budget
    print('WFC/ADC reflec. (WS)', np.interp(lam1, pfue_ws['Wavelength'] * 10000, pfue_ws['Reflexion']))
    print('WFC/ADC trans. (WS)', np.interp(lam1, pfue_ws['Wavelength'] * 10000, pfue_ws['Transmission']))
    print('WFC/ADC vign. (WS)', np.interp(lam1, pfue_ws['Wavelength'] * 10000, pfue_ws['Vignetting']))
    print('WFC/ADC total (WS)', np.interp(lam1, pfue_ws['Wavelength'] * 10000, pfue_ws['Total']))
    pfue_thr_ws = np.interp(lam, pfue_ws['Wavelength'] * 10000, pfue_ws['Total'])

    # plot and save
    plt.figure(figsize=(10, 5))
    # add Will's numbers to the plot
    plt.semilogx(pfue_ws['Wavelength'] * 1000., pfue_ws['Reflexion'], 'b', label="Reflection")
    plt.plot(pfue_ws['Wavelength'] * 1000., pfue_ws['Transmission'], 'r', label="Transmission")
    plt.plot(pfue_ws['Wavelength'] * 1000., pfue_ws['Vignetting'], 'g', label="Vignetting")
    plt.plot(pfue_ws['Wavelength'] * 1000., pfue_ws['Total'], 'k', label="Total")
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_pfue.png', bbox_inches="tight")

    # = = = = = = = = =
    # Positioners (FRD)
    # = = = = = = = = =

    # FRD losses
    poss_frd = np.ones_like(lam) * 0.97
    # save file
    poss_data = Table([lam, poss_frd], names=['lamA', 'thr'])
    ascii.write(poss_data, path + 'mse_itc_throughput_poss.dat', overwrite=True)
    # print values at Lam0 and Lam1 for budget
    print('PosS (FRD)', np.interp(lam1, lam, poss_frd))

    # = = = = = = = = =
    # Fibre train
    # = = = = = = = = =

    #  FRD losses
    fib_frd = np.ones_like(lam) * 0.95

    #  Fresnel losses
    fresn_hr = np.ones_like(lam) * 0.97
    fresn_hr[lam < 4000] = 0.96
    fresn_hr[lam > 9000] = 0.
    fresn_lmr = np.ones_like(lam) * 0.88
    fresn_lmr[lam < 10000] = 0.94

    # Transmission
    # Fibre transmission (db/km)
    polyfbp = Table.read(path + 'polymicro_fbp.csv', format='csv')
    # length (in km)
    lmr_len = 30 / 1e3
    hr_len = 50 / 1e3
    lmr_len_b = 50 / 1e3
    hr_len_b = 30 / 1e3
    # attenuation as a power ratio(see definition of dB)
    trans_lmr = 10. ** (- polyfbp['att'].data * lmr_len / 10.)
    trans_hr = 10. ** (- polyfbp['att'].data * hr_len / 10.)
    trans_lmr_b = 10. ** (- polyfbp['att'].data * lmr_len_b / 10.)
    trans_hr_b = 10. ** (- polyfbp['att'].data * hr_len_b / 10.)
    # interpol
    trans_lmr = np.interp(lam, polyfbp['wave'].data * 10, trans_lmr)
    trans_hr = np.interp(lam, polyfbp['wave'].data * 10, trans_hr)
    trans_lmr_b = np.interp(lam, polyfbp['wave'].data * 10, trans_lmr_b)
    trans_hr_b = np.interp(lam, polyfbp['wave'].data * 10, trans_hr_b)

    # Print values at reference wavelengths
    print('FiTS LMR - Trans FBP', np.interp(lam1, lam, trans_lmr))
    print('FiTS HR - Trans FBP', np.interp(lam1, lam, trans_hr))

    # Fibre transmission (db/km) - FBPI
    polyfbpi = Table.read(path + 'polymicro_fbpi.csv', format='csv')
    # y-axis had different scales on the figure so we combine them here
    polyfbpi['att'].data[polyfbpi['wave'] < 397.5] *= 1000  # was in db/m, scale 0 to 2
    polyfbpi['att'].data[polyfbpi['wave'] > 397.5] *= 25  # was in db/km, scale 0 to 50
    # attenuation as a power ratio(see definition of dB)
    att_lmr_fbpi = 10. ** (- polyfbpi['att'].data * lmr_len / 10.)
    att_hr_fbpi = 10. ** (- polyfbpi['att'].data * hr_len / 10.)
    # interpol
    trans_lmr_fbpi = np.interp(lam, polyfbpi['wave'].data * 10, att_lmr_fbpi)
    trans_hr_fbpi = np.interp(lam, polyfbpi['wave'].data * 10, att_hr_fbpi)

    plt.figure(figsize=(10, 5))
    plt.semilogx(lam / 10., trans_lmr_fbpi, 'b:', label='LMR - FBPI (30 m)')
    plt.plot(lam / 10., trans_hr_fbpi, 'r:', label='HR - FBPI (50 m)')

    # Print values at reference wavelengths
    print('FiTS LMR - Trans FBPI', np.interp(lam1, lam, trans_lmr_fbpi))
    print('FiTS HR - Trans FBPI', np.interp(lam1, lam, trans_hr_fbpi))

    # Total fiber
    fits_lmr_thr = trans_lmr * fresn_lmr * fib_frd
    fits_hr_thr = trans_hr * fresn_hr * fib_frd
    fits_lmr_thr_b = trans_lmr_b * fresn_lmr * fib_frd
    fits_hr_thr_b = trans_hr_b * fresn_hr * fib_frd
    fits_lmr_thr_fbpi = trans_lmr_fbpi * fresn_lmr * fib_frd
    fits_hr_thr_fbpi = trans_hr_fbpi * fresn_hr * fib_frd
    # save file
    fits_lmr_data = Table([lam, fits_lmr_thr], names=['lamA', 'thr'])
    ascii.write(fits_lmr_data, path + 'mse_itc_throughput_fits_lmr.dat', overwrite=True)
    fits_hr_data = Table([lam, fits_hr_thr], names=['lamA', 'thr'])
    ascii.write(fits_hr_data, path + 'mse_itc_throughput_fits_hr.dat', overwrite=True)
    # Print values at reference wavelengths
    print('FiTS LMR', np.interp(lam1, lam, fits_lmr_thr))
    print('FiTS HR', np.interp(lam1, lam, fits_hr_thr))
    print('FiTS LMR - FBPI', np.interp(lam1, lam, fits_lmr_thr_fbpi))
    print('FiTS HR - FBPI', np.interp(lam1, lam, fits_hr_thr_fbpi))

    # plot and save
    plt.figure(figsize=(10, 5))
    plt.semilogx(lam / 10., fits_lmr_thr, 'b', label='LMR - FBP (30 m)')
    plt.plot(lam / 10., fits_hr_thr, 'r', label='HR - FBP (50 m)')
    plt.semilogx(lam / 10., fits_lmr_thr_b, 'b--', label='LMR - FBP (50 m)')
    plt.plot(lam / 10., fits_hr_thr_b, 'r--', label='HR - FBP (30 m)')
    plt.semilogx(lam / 10., fits_lmr_thr_fbpi, 'b:', label='LMR - FBPI (30 m)')
    plt.plot(lam / 10., fits_hr_thr_fbpi, 'r:', label='HR - FBPI (50 m)')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_fits.png', bbox_inches="tight")

    # = = = = = = = = =
    # Spectrographs
    # = = = = = = = = =

    margin_factor = 0.8  # 1.0 to use CoDR numbers directly, <1.0 to add a margin

    # --- NEW SPECTROGRAPH(7 / 27 / 2017) - 20 % below curve from LMR / HR CoDR
    # LR
    spec_lr_1 = Table.read(path + 'lr_throughput_w1.csv', format='csv')
    spec_lr_2 = Table.read(path + 'lr_throughput_w2.csv', format='csv')
    spec_lr_3 = Table.read(path + 'lr_throughput_w3.csv', format='csv')
    spec_lr_4 = Table.read(path + 'lr_throughput_w4.csv', format='csv')
    spec_lr_5 = Table.read(path + 'lr_throughput_w5.csv', format='csv')
    spec_lr = np.empty_like(lam).reshape(1500, 1).dot(np.zeros((1, 5)))
    spec_lr[(lam <= 5600) & (lam >= 3600), 0] = np.interp(lam[(lam <= 5600) & (lam >= 3600)],
                                                          spec_lr_1['wave'] * 10, spec_lr_1['thr_arm1']) * margin_factor
    spec_lr[(lam <= 7400) & (lam >= 5400), 1] = np.interp(lam[(lam <= 7400) & (lam >= 5400)],
                                                          spec_lr_2['wave'] * 10, spec_lr_2['thr_arm2']) * margin_factor
    spec_lr[(lam <= 9850) & (lam >= 7150), 2] = np.interp(lam[(lam <= 9850) & (lam >= 7150)],
                                                          spec_lr_3['wave'] * 10, spec_lr_3['thr_arm3']) * margin_factor
    spec_lr[(lam <= 13200) & (lam >= 9600), 3] = np.interp(lam[(lam <= 13200) & (lam >= 9600)],
                                                           spec_lr_4['wave'] * 10, spec_lr_4['thr_arm4']) * margin_factor
    spec_lr[(lam <= 18000) & (lam >= 14570), 4] = np.interp(lam[(lam <= 18000) & (lam >= 14570)],
                                                            spec_lr_5['wave'] * 10, spec_lr_5['thr_arm5']) * margin_factor

    # MR
    spec_mr_1 = Table.read(path + 'mr_throughput_w1.csv', format='csv')
    spec_mr_2 = Table.read(path + 'mr_throughput_w2.csv', format='csv')
    spec_mr_3 = Table.read(path + 'mr_throughput_w3.csv', format='csv')
    spec_mr_nograt_1 = Table.read(path + 'mr_throughput_nograt_w1.csv', format='csv')
    spec_mr_nograt_2 = Table.read(path + 'mr_throughput_nograt_w2.csv', format='csv')
    spec_mr_nograt_3 = Table.read(path + 'mr_throughput_nograt_w3.csv', format='csv')
    # without grating (with 20% margin)
    spec_mr_nograt = np.empty_like(lam).reshape(1500, 1).dot(np.zeros((1, 3)))
    #spec_mr_nograt[(lam <= 5100) & (lam >= 3910), 0] = np.interp(lam[(lam <= 5100) & (lam >= 3910)], spec_mr_nograt_1['wav'] * 10, spec_mr_nograt_1['thr'] * 0.8)
    #spec_mr_nograt[(lam <= 7000) & (lam >= 5760), 1] = np.interp(lam[(lam <= 7000) & (lam >= 5760)], spec_mr_nograt_2['wav'] * 10, spec_mr_nograt_2['thr'] * 0.8)
    #spec_mr_nograt[(lam <= 9000) & (lam >= 7370), 2] = np.interp(lam[(lam <= 9000) & (lam >= 7370)], spec_mr_nograt_3['wav'] * 10, spec_mr_nograt_3['thr'] * 0.8)

    spec_mr_nograt[:, 0] = np.interp(lam, spec_mr_nograt_1['wav'] * 10, spec_mr_nograt_1['thr'] * margin_factor)
    spec_mr_nograt[:, 1] = np.interp(lam, spec_mr_nograt_2['wav'] * 10, spec_mr_nograt_2['thr'] * margin_factor)
    spec_mr_nograt[:, 2] = np.interp(lam, spec_mr_nograt_3['wav'] * 10, spec_mr_nograt_3['thr'] * margin_factor)

    # full throughput first (with margin)
    spec_mr = np.empty_like(lam).reshape(1500, 1).dot(np.zeros((1, 4)))
    spec_mr[(lam <= 5100) & (lam >= 3910), 0] = np.interp(lam[(lam <= 5100) & (lam >= 3910)],
                                                          spec_mr_1['wave'] * 10, spec_mr_1['thr_arm1']) * margin_factor
    spec_mr[(lam <= 7000) & (lam >= 5760), 1] = np.interp(lam[(lam <= 7000) & (lam >= 5760)],
                                                          spec_mr_2['wave'] * 10, spec_mr_2['thr_arm2']) * margin_factor
    spec_mr[(lam <= 9000) & (lam >= 7370), 2] = np.interp(lam[(lam <= 9000) & (lam >= 7370)],
                                                          spec_mr_3['wave'] * 10, spec_mr_3['thr_arm3']) * margin_factor
    spec_mr[(lam <= 18000) & (lam >= 14570), 3] = np.interp(lam[(lam <= 18000) & (lam >= 14570)],
                                                            spec_lr_5['wave'] * 10, spec_lr_5['thr_arm5']) * margin_factor

    # HR
    spec_hr_nograt = np.empty_like(lam).reshape(1500, 1).dot(np.zeros((1, 3)))
    # without gratings first
    # blue
    spec_hr_nograt_blue = Table.read(path + 'hr_throughput_nogratings_blue.csv', format='csv')
    sel = (lam <= np.max(spec_hr_nograt_blue['wav']) * 10) & (lam >= np.min(spec_hr_nograt_blue['wav']) * 10)
    spec_hr_nograt[sel, 0] = np.interp(lam[sel], spec_hr_nograt_blue['wav'] * 10, spec_hr_nograt_blue['throughput'] * margin_factor)
    # green
    spec_hr_nograt_green = Table.read(path + 'hr_throughput_nogratings_green.csv', format='csv')
    sel = (lam <= np.max(spec_hr_nograt_green['wav']) * 10) & (lam >= np.min(spec_hr_nograt_green['wav']) * 10)
    spec_hr_nograt[sel, 1] = np.interp(lam[sel], spec_hr_nograt_green['wav'] * 10, spec_hr_nograt_green['throughput'] * margin_factor)
    # red
    spec_hr_nograt_red = Table.read(path + 'hr_throughput_nogratings_red.csv', format='csv')
    sel = (lam <= np.max(spec_hr_nograt_red['wav']) * 10) & (lam >= np.min(spec_hr_nograt_red['wav']) * 10)
    spec_hr_nograt[sel, 2] = np.interp(lam[sel], spec_hr_nograt_red['wav'] * 10, spec_hr_nograt_red['throughput'] * margin_factor)
    # gratings alone
    spec_hr_grat = np.empty_like(lam).reshape(1500, 1).dot(np.zeros((1, 3)))
    # blue
    spec_hr_grat_blue = Table.read(path + 'HR_throughput_gratings_blue.csv', format='csv')
    sel = (lam <= np.max(spec_hr_grat_blue['wav']) * 10) & (lam >= np.min(spec_hr_grat_blue['wav']) * 10)
    spec_hr_grat[sel, 0] = np.interp(lam[sel], spec_hr_grat_blue['wav'] * 10, spec_hr_grat_blue['throughput'] / 100)
    # green
    spec_hr_grat_green = Table.read(path + 'HR_throughput_gratings_green.csv', format='csv')
    sel = (lam <= np.max(spec_hr_grat_green['wav']) * 10) & (lam >= np.min(spec_hr_grat_green['wav']) * 10)
    spec_hr_grat[sel, 1] = np.interp(lam[sel], spec_hr_grat_green['wav'] * 10, spec_hr_grat_green['throughput'] / 100)
    # red
    spec_hr_grat_red = Table.read(path + 'HR_throughput_gratings_red.csv', format='csv')
    sel = (lam <= np.max(spec_hr_grat_red['wav']) * 10) & (lam >= np.min(spec_hr_grat_red['wav']) * 10)
    spec_hr_grat[sel, 2] = np.interp(lam[sel], spec_hr_grat_red['wav'] * 10, spec_hr_grat_red['throughput'] / 100)
    # combine
    spec_hr = spec_hr_nograt * spec_hr_grat

    # save file
    spec_lr_data = Table([lam, spec_lr[:, 0], spec_lr[:, 1], spec_lr[:, 2], spec_lr[:, 3], spec_lr[:, 4]],
                         names=['lamA', 'thr_arm1', 'thr_arm2', 'thr_arm3', 'thr_arm4', 'thr_arm5'])
    ascii.write(spec_lr_data, path + 'mse_itc_throughput_spec_lr.dat', overwrite=True)
    spec_mr_data = Table([lam, spec_mr[:, 0], spec_mr[:, 1], spec_mr[:, 2], spec_mr[:, 3]],
                         names=['lamA', 'thr_arm1', 'thr_arm2', 'thr_arm3', 'thr_arm4'])
    ascii.write(spec_mr_data, path + 'mse_itc_throughput_spec_mr.dat', overwrite=True)
    spec_hr_data = Table([lam, spec_hr[:, 0], spec_hr[:, 1], spec_hr[:, 2]],
                         names=['lamA', 'thr_arm1', 'thr_arm2', 'thr_arm3'])
    ascii.write(spec_hr_data, path + 'mse_itc_throughput_spec_hr.dat', overwrite=True)

    # Print values at reference wavelengths
    for i in range(5):
        print('Spectro LR', np.interp(lam1, lam, spec_lr[:, i]))
    for i in range(3):
        print('Spectro MR (no grat)', np.interp(lam1, lam, spec_mr_nograt[:, i]))
    for i in range(3):
        print('Spectro MR', np.interp(lam1, lam, spec_mr[:, i]))
    for i in range(3):
        print('Spectro HR (no grat)', np.interp(lam1, lam, spec_hr_nograt[:, i]))
    for i in range(3):
        print('Spectro HR ', np.interp(lam1, lam, spec_hr[:, i]))

    # plot and save - LR
    plt.figure(figsize=(10, 5))
    for i in range(5):
        plt.semilogx(lam / 10., spec_lr[:, i], label='LR '+['B','G','R','J','H'][i], c=['b','g','r','#AA0000','#660000'][i])
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_specLR.png', bbox_inches="tight")
    # plot and save - MR
    plt.figure(figsize=(10, 5))
    for i in range(3):
        plt.semilogx(lam / 10., spec_mr[:, i], label='MR ' + ['B', 'G', 'R'][i], c=['b', 'g', 'r'][i], linestyle='-')
    for i in range(3):
        plt.semilogx(lam / 10., spec_mr_nograt[:, i], label='MR ' + ['B', 'G', 'R'][i] + ' (no grat.)',
                     c=['b', 'g', 'r'][i], linestyle='--')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_specMR.png', bbox_inches="tight")
    # plot and save - HR
    plt.figure(figsize=(10, 5))
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr[:, i], label='HR '+['B','G','R'][i], c=['b','g','r'][i], linestyle='-')
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr_nograt[:, i], label='HR '+['B','G','R'][i]+' (no grat.)', c=['b','g','r'][i], linestyle='--')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_specHR.png', bbox_inches="tight")
    # plot and save - HR (gratings only)
    plt.figure(figsize=(10, 5))
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr[:, i], label='HR ' + ['B', 'G', 'R'][i], c=['b', 'g', 'r'][i], linestyle='-')
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr_nograt[:, i], label='HR ' + ['B', 'G', 'R'][i] + ' (no grat.)',
                     c=['b', 'g', 'r'][i], linestyle='--')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_specHR.png', bbox_inches="tight")

    # = = = = = = = = =
    # Overall throughput
    # = = = = = = = = =

    # Print values at reference wavelengths
    for i in range(5):
        print('MSE LR', np.interp(lam1, lam, spec_lr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr))
    for i in range(4):
        print('MSE MR', np.interp(lam1, lam, spec_mr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr))
    for i in range(3):
        print('MSE MR (no grat)', np.interp(lam1, lam, spec_mr_nograt[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr))
    for i in range(3):
        print('MSE HR', np.interp(lam1, lam, spec_hr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_hr_thr))
    for i in range(3):
        print('MSE HR (no grat)',
              np.interp(lam1, lam, spec_hr_nograt[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_hr_thr))

    # plot and save - LR
    plt.figure(figsize=(10, 5))
    for i in range(5):
        plt.semilogx(lam / 10., spec_lr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr, label='LR '+['B','G','R','J','H'][i], c=['b','g','r','#AA0000','#660000'][i])
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_LR.png', bbox_inches="tight")
    # plot and save - MR
    plt.figure(figsize=(10, 5))
    for i in range(3):
        plt.semilogx(lam / 10., spec_mr_nograt[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr,
                     label='MR ' + ['B', 'G', 'R'][i] + ' (no grat.)', c=['b', 'g', 'r'][i], linestyle='--')
    for i in range(3):
        plt.semilogx(lam / 10., spec_mr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_lmr_thr,
                     label='MR ' + ['B', 'G', 'R'][i], c=['b', 'g', 'r'][i], linestyle='-')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_MR.png', bbox_inches="tight")
    # plot and save - HR
    plt.figure(figsize=(10, 5))
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_hr_thr, label='HR '+['B','G','R'][i], c=['b','g','r'][i], linestyle='-')
    for i in range(3):
        plt.semilogx(lam / 10., spec_hr_nograt[:, i] * struc_thr * m1_thr * pfue_thr_ws * poss_frd * fits_hr_thr, label='HR '+['B','G','R'][i]+' (no grat.)', c=['b','g','r'][i], linestyle='--')
    plt.plot([360, 360, 370, 370, 400, 400, 1800, 1800], [0, 10, 10, 0, 0, 10, 10, 0], 'c:')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.xlim([360, 1800])
    plt.ylim([0.0, 1.05])
    plt.legend()
    plt.savefig(path + 'mse_itc_throughput_HR.png', bbox_inches="tight")

mse_itc_throughput()
print('Done!')