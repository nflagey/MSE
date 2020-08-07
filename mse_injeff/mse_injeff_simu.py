#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
import os
import fnmatch
from astropy.io import fits
from matplotlib import pyplot as plt
import scipy.io
import photutils as phot
import time
from scipy import signal

# Some default parameters
plate_scale = 106.7e-6  # meter per arcsec
pix_scale = 0.3e-6  # meter per pixel (Zemax)
fov_rad_deg = 0.76  # degrees (from optical design TN)
fov_rad_m = plate_scale * fov_rad_deg * 3600
spine_length = 300e-3  # meter

# Foci
foci = {'N100': -100e-6, 'N25': -25e-6, '0': 0, 'P25': 25e-6, 'P100': 100e-6}
nfoci = len(foci)
#  Wavelengths
wave = np.asarray([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1784])  # in nm
nwave = len(wave)
# Weights for optimization
wave_weight = {'LR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 4, 2, 4, 2, 4, 1]),
               'MR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 3, 2, 1, 0, 0, 0]),
               'HR': np.interp(wave, [360, 370, 550, 750, 1000, 1300, 1800], [1, 3, 2, 1, 0, 0, 0])}

# Injection efficiency budget
# Number of simulations in each configuration
nsimu = 50
# Number of steps in each observation
nstep = 25
# Positioners parfocality
z_parf_poss = 25e-6           # 25um max - good with SysCoDR budget
# Fiber in positioners parfocality (addition)
z_parf_fits = 5e-6            # 5um max - good with SysCoDR budget
# AIV misalignment
z_aiv = 50e-6                 # 50um max - good with SysCoDR budget
# Residual errors in lookup tables
z_lookup = 10e-6        # 10um max - good with SysCoDR budget
# Gravity loads variations during exposure
z_grav = 1e-6                 # 1um max
# Thermal loads variations during exposure
r_therm = 14e-6               # 14um max at edge of field (i.e. +/- 7 um)
z_therm = 2e-6                # 2um max
# InRo tilt
z_inro = 30e-6                # 30um max at edge of field
# Hexapod correction accuracy
z_pfhs = 5e-6                 # 5um max
# Target coordinates error
r_tgtcoo = 10e-6              # 10um max
# Sky to focal surface mapping
r_skymap = 5e-6               # 5um max - good with SysCoDR budget
# PosS accuracy
r_poss_budget = 6e-6          # 6um rms - good with SysCoDR budget
r_poss_aao = [4.5e-6, 8.1e-6]   # rms, max
r_poss_ustc = [5.1e-6, 21.7e-6]  # rms, max
# Guiding error
r_guid = 10e-6              # 10um rms (Kevin Ho 5/31/2017) ... but 2um rms in IE budget
# Instrument Rotator accuracy --> gaussian spread (on an arc)
r_inro = 5e-6               # 5um rms at edge of field
# Plate scale variations
r_plate = 10e-6             # 10um max at edge of field
# ADR residuals
r_adr = 15e-6               # 30um max at edge of field (should depend on ZD) but aiming for center so +/-15um
# Vibrations
r_vibr = 1e-6               # 1um max


class MsePSF:

    def __init__(self, design='segments', src_type='point', zenith='00', field='X_+0.530', seeing=0.4, spectro='LR', fibdiam=1.0):

        # Design ('segments' or 'no_segments')
        self.design = design
        # Target
        self.src_type = src_type
        # Observing conditions
        self.seeing = seeing  # at 500 nm
        self.spectro = spectro
        self.fibdiam = fibdiam
        # Zemax parameters
        self.field = field
        self.zenith = zenith

        # Cube and header
        self.cube = None
        self.header = None
        # Cube convolved with IQ
        self.cubeconvfile = 'results/IQ_cubes/' + design + '_IQ_cube_Z_' + self.zenith + '_Field_' + self.field + '_IQ_' + str(seeing) + '.fits'
        # Optimal position and optimal IE curve
        self.optim = None
        self.optim_ie = None
        # Offsets
        self.x_off = None
        self.y_off = None
        self.z_off = None
        # IE curve
        self.ie_simu = None
        self.ie_simu_stddev = None

    def open_zemax(self):
        """ Open Zemax file (segments combined, one wavelength per frame, all focus values)
        """

        zemax_data = []
        zemax_header = []
        for focus in foci:
            path = '/Users/nflagey/WORK/MSE/MSE_Injection/PSFs/' + self.design + '/3_cubes/' \
                   + 'Z_' + self.zenith + '/Focus_' + focus + '/Field_' + self.field
            cube_file = os.listdir(path)
            hdu = fits.open(os.path.join(path, cube_file[0]))
            zemax_data.append(hdu[0].data)
            zemax_header.append(hdu[0].header)
            self.cube = zemax_data
            self.header = zemax_header

    def apply_iq(self):
        """ Apply natural seeing and others to Zemax PSF
        """

        if self.seeing == 0:
            # Return the same (just change format)
            cubeconv = np.zeros_like(self.cube)
            for w in range(nwave):
                for f in range(nfoci):
                    cubeconv[f, w, :, :] = self.cube[f][w, :, :]

        else:
            # Seeing (fwhm in arcsec) as a function of wavelength
            seeing_wave = (wave / 500.) ** (-1. / 5) * self.seeing

            # Output
            cubeconv = np.zeros_like(self.cube)

            # Loop over wavelength frames to build kernel
            for w in range(nwave):
                kernel = ie_moffat(seeing_wave[w])
                # Convolution (loop over foci)
                for f in range(nfoci):
                    cubeconv[f, w, :, :] = signal.fftconvolve(self.cube[f][w, :, :], kernel, mode='same')
            # Save into a file
            hdu = fits.PrimaryHDU(cubeconv)
            hdu.writeto(self.cubeconvfile, overwrite=True)

    def find_optim(self):
        """ Find optimal position on the zemax data after IQ
        """

        # Open FITS file
        hdu = fits.open(self.cubeconvfile)
        cubeconv = hdu[0].data
        hdr = hdu[0].header

        # In-focus only
        f0 = [i for i in range(nfoci) if list(foci.values())[i] == 0][0]

        # Find peak of emission at each wavelength
        peaks = np.zeros((2, nwave), dtype=int)
        for k in range(nwave):
            peak = np.where(cubeconv[f0, k, :, :] == np.max(cubeconv[f0, k, :, :]))
            peaks[:, k] = peak[0]

        # Loop in a scarcely sample box near peaks
        skip = 5
        nskip = 3
        x = range(np.min(peaks[0, :]) - nskip * skip, np.max(peaks[0, :]) + nskip * skip + 1, skip)
        y = range(np.min(peaks[1, :]) - nskip * skip, np.max(peaks[1, :]) + nskip * skip + 1, skip)

        # Create all the apertures
        aper = phot.CircularAperture([(xi, yj) for xi in x for yj in y], r=self.fibdiam / 2 * plate_scale / pix_scale)

        # Aperture photometry on all apertures and all wavelengths
        res = [phot.aperture_photometry(cubeconv[f0, k, :, :] / np.sum(self.cube[f0][k, :, :]), aper)
               for k in range(nwave)]
        # This gives 11 "tables" of aperture photometry object, each table being nx*ny values

        # Weighted average by hand because we are not dealing with a simple np.array
        res_weighted = [np.sum([res[i]['aperture_sum'].quantity[j] * wave_weight[self.spectro][i]
                                for i in range(nwave)])/np.sum(wave_weight[self.spectro])
                        for j in range(len(x)*len(y))]
        # This gives a nx*ny list of values of weighted averaged IE

        # Find maximum value
        optim = np.where(res_weighted == np.max(res_weighted))
        # Find position of maximum
        optim_pos = [int(res[0]['xcenter'].value[optim]), int(res[0]['ycenter'].value[optim])]
        # Now we have the optimal position over a coarse grid

        # Plot
        plt.figure()
        plt.imshow(np.resize(res_weighted, (len(x), len(y))))
        plt.colorbar()
        plt.savefig('results/optim/' + self.design + '_optim_coarse_' + self.spectro + '_fib' + str(self.fibdiam) +
                    '_Z' + self.zenith + '_Field_' + self.field + '_Seeing_' + str(self.seeing) + '.png')
        plt.close()

        # Loop in a thinly sample but smaller box near peaks
        x = range(optim_pos[0] - skip, optim_pos[0] + skip)
        y = range(optim_pos[1] - skip, optim_pos[1] + skip)

        # Create all the apertures
        aper = phot.CircularAperture([(xi, yj) for xi in x for yj in y], r=self.fibdiam / 2 * plate_scale / pix_scale)

        # Aperture photometry on all apertures and all wavelengths
        res = [phot.aperture_photometry(cubeconv[f0, k, :, :] / np.sum(self.cube[f0][k, :, :]), aper)
               for k in range(nwave)]

        # This gives 11 "tables" of aperture photometry object, each table being nx*ny values

        # Weighted average by hand because we are not dealing with a simple np.array
        res_weighted = [np.sum([res[i]['aperture_sum'].quantity[j] * wave_weight[self.spectro][i]
                                for i in range(nwave)])/np.sum(wave_weight[self.spectro])
                        for j in range(len(x)*len(y))]
        # This gives a nx*ny list of values of weighted averaged IE

        # Find maximum value
        optim = np.where(res_weighted == np.max(res_weighted))
        # Find position of maximum
        optim_pos = [int(res[0]['xcenter'].value[optim]), int(res[0]['ycenter'].value[optim])]
        # Now we have the optimal position over a fine grid

        # Plot
        plt.figure()
        plt.imshow(np.resize(res_weighted, (len(x), len(y))))
        plt.colorbar()
        plt.savefig('results/optim/' + self.design + '_optim_fine_' + self.spectro + '_fib' + str(self.fibdiam) +
                    '_Z' + self.zenith + '_Field_' + self.field + '_Seeing_' + str(self.seeing) + '.png')
        plt.close()

        # Save optimal position and optimal IE into dictionary
        self.optim = optim_pos
        self.optim_ie = np.asarray([res[i]['aperture_sum'].quantity[optim].value for i in range(nwave)])

        # Save optimal position into header of FITS file so IE can find it
        hdr['OPTX' + self.spectro[0] + str(self.fibdiam)[0] + str(self.fibdiam)[2:4]] = optim_pos[0]
        hdr['OPTY' + self.spectro[0] + str(self.fibdiam)[0] + str(self.fibdiam)[2:4]] = optim_pos[1]
        hdu.writeto(self.cubeconvfile, overwrite=True)

    def ie_budget(self):
        """ Pick up random values for XYZ offsets according to IE budget
        """
        # ===== Position error that is fixed for a given observation or given positioner

        # As delivered
        #      parfocality of positioners tips - MAX
        z_parf_poss_samp = np.random.random_sample(nsimu) * 2 * z_parf_poss - z_parf_poss
        #      additional scatter due to FiTS - MAX
        z_parf_fits_samp = np.random.random_sample(nsimu) * 2 * z_parf_fits - z_parf_fits

        # AIV misalignment - MAX
        z_aiv_samp = np.random.random_sample(nsimu) * 2 * z_aiv - z_aiv

        # Deformation of structure, lookup table residuals
        #   Gravity and thermal loads - MAX
        z_lookup_samp = np.random.random_sample(nsimu) * 2 * z_lookup - z_lookup

        # Relative fiber positioning
        #   Target location - MAX distance, random direction
        r_tgtcoo_samp = np.random.random_sample(nsimu) * r_tgtcoo
        t_tgtcoo_samp = np.random.random_sample(nsimu) * 2. * np.pi
        x_tgtcoo_samp = r_tgtcoo_samp * np.cos(t_tgtcoo_samp)
        y_tgtcoo_samp = r_tgtcoo_samp * np.sin(t_tgtcoo_samp)
        #   Sky to focal surface mapping - MAX distance, random direction
        r_skymap_samp = np.random.random_sample(nsimu) * r_skymap
        t_skymap_samp = np.random.random_sample(nsimu) * 2. * np.pi
        x_skymap_samp = r_skymap_samp * np.cos(t_skymap_samp)
        y_skymap_samp = r_skymap_samp * np.sin(t_skymap_samp)
        #   Positioner accuracy - RMS
        r_poss_samp = np.random.normal(size=nsimu) * r_poss_aao[0]
        r_poss_samp[np.abs(r_poss_samp) > r_poss_aao[1]] = r_poss_aao[1]
        t_poss_samp = np.random.random_sample(nsimu) * 2 * np.pi
        x_poss_samp = r_poss_samp * np.cos(t_poss_samp)
        y_poss_samp = r_poss_samp * np.sin(t_poss_samp)

        # Fiber out of focus position
        #   Fiber tilt - need distribution
        z_defoc_poss_samp = self.simu_defoc_poss()

        # ===== Position error that varies during a given observation

        # InRo tilt - 30 microns MAX, this is varying slowly and
        #                  following a predefined path so we can pre-compute
        #                  all values
        z_inro_samp = np.transpose([self.simu_defoc_inrotilt() for j in range(nsimu)])

        # Gravity loads (Z_ - MAX
        z_grav_samp = [np.random.random_sample(nstep) * z_grav * 2 - z_grav for j in range(nsimu)]
        flip = [np.round(np.random.random_sample(1)) for j in range(nsimu)]  # flip a coin to increase/decrease
        z_grav_samp = np.sort(z_grav_samp)
        z_grav_samp = np.transpose([z_grav_samp[j, ::-1] if flip[j] else np.sort(z_grav_samp[j, :]) for j in range(nsimu)])

        # Thermal loads (Z) - MAX
        z_therm_samp = [np.random.random_sample(nstep) * z_therm * 2 - z_therm for j in range(nsimu)]
        flip = [np.round(np.random.random_sample(1)) for j in range(nsimu)]  # flip a coin to increase/decrease
        z_therm_samp = np.sort(z_therm_samp)
        z_therm_samp = np.transpose([z_therm_samp[j, ::-1] if flip[j] else np.sort(z_therm_samp[j, :]) for j in range(nsimu)])

        # Hexapod correction residual (Z) - MAX
        z_pfhs_samp = np.transpose([np.random.random_sample(nstep) * z_pfhs * 2 - z_pfhs for j in range(nsimu)])

        # Vibrations (XY) - MAX
        r_vibr_samp = np.transpose([np.random.random_sample(nstep) * r_vibr for j in range(nsimu)])
        t_vibr_samp = np.transpose([np.random.random_sample(nstep) * 2 * np.pi for j in range(nsimu)])
        x_vibr_samp = r_vibr_samp * np.cos(t_vibr_samp)
        y_vibr_samp = r_vibr_samp * np.sin(t_vibr_samp)

        # Thermal loads (XY) - MAX
        r_therm_samp = [np.random.random_sample(nstep) * r_therm * 2 - r_therm for j in range(nsimu)]
        theta_therm_samp = [np.random.random_sample(nstep) * 2 * np.pi for j in range(nsimu)]
        x_therm_samp = r_therm_samp * np.cos(theta_therm_samp)
        y_therm_samp = r_therm_samp * np.sin(theta_therm_samp)
        flip = [np.round(np.random.random_sample(1)) for j in range(nsimu)]  # flip a coin to increase/decrease
        x_therm_samp = np.sort(x_therm_samp)
        x_therm_samp = np.transpose([x_therm_samp[j, ::-1] if flip[j] else np.sort(x_therm_samp[j, :])
                                     for j in range(nsimu)])
        flip = [np.round(np.random.random_sample(1)) for j in range(nsimu)]  # flip a coin to increase/decrease
        y_therm_samp = np.sort(y_therm_samp)
        y_therm_samp = np.transpose([y_therm_samp[j, ::-1] if flip[j] else np.sort(y_therm_samp[j, :])
                                     for j in range(nsimu)])

        # Telescope motion during exposure (XY) - RMS
        #      Guiding errors
        r_guid_samp = np.transpose([np.random.normal(size=nstep) * r_guid for j in range(nsimu)])
        theta_guid_samp = np.transpose([np.random.random_sample(nstep) * 2 * np.pi for j in range(nsimu)])
        x_guid_samp = r_guid_samp * np.cos(theta_guid_samp)
        y_guid_samp = r_guid_samp * np.sin(theta_guid_samp)
        #      Rotation rate - 5 microns RMS at edge of field, so does depend
        #                      on target position in FoV
        inro_samp = np.transpose([self.simu_offset_inro() for j in range(nsimu)])
        x_inro_samp, y_inro_samp = inro_samp[:, 0, :], inro_samp[:, 1, :]

        # Atmospheric differential refraction (XY) - special
        #      WFC/ADC distortion + ADR drift
        #      - or - Open-loop accuracy
        # x_adr_samp, y_adr_samp = np.transpose([self.simu_offset_adr() for j in range(nsimu)])
        adr_samp = np.transpose([self.simu_offset_adr() for j in range(nsimu)])
        x_adr_samp, y_adr_samp = adr_samp[:, 0, :], adr_samp[:, 1, :]

        # Plate scale variations (XY) - special
        #      thermal and gravity variations - this should be radial in the
        #                                       whole FoV, not radial for each fiber
        plate_samp = np.transpose([self.simu_offset_plate() for j in range(nsimu)])
        x_plate_samp, y_plate_samp = plate_samp[:, 0, :], plate_samp[:, 1, :]

        # Sum all offsets properly (some NSIMU, some NSTEP)
        x_off = np.asarray([np.asarray([x_tgtcoo_samp[i] + x_skymap_samp[i] + x_poss_samp[i] for i in range(nsimu)])
                            for j in range(nstep)] + x_vibr_samp + x_guid_samp + x_adr_samp + x_plate_samp
                           + x_inro_samp + x_therm_samp)

        y_off = np.asarray([np.asarray([y_tgtcoo_samp[i] + y_skymap_samp[i] + y_poss_samp[i] for i in range(nsimu)])
                            for j in range(nstep)] + y_vibr_samp + y_guid_samp + y_adr_samp + y_plate_samp
                           + y_inro_samp + y_therm_samp)

        z_off = np.asarray([np.asarray([z_parf_poss_samp[i] + z_parf_fits_samp[i] + z_aiv_samp[i] + z_lookup_samp[i]
                                        + z_defoc_poss_samp[i] for i in range(nsimu)]) for j in range(nstep)]
                           + z_inro_samp + z_pfhs_samp + z_therm_samp + z_grav_samp)

        # === PLOTS ===
        # --- Plot offsets ---
        plt.figure(figsize=(10, 8))
        plt.rc('font', size=6)  # controls default text sizes
        plt.rc('axes', titlesize=6)  # fontsize of the axes title
        plt.rc('axes', labelsize=6)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=6)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=6)  # fontsize of the tick labels
        plt.rc('legend', fontsize=6)  # legend fontsize
        plt.rc('figure', titlesize=8)  # fontsize of the figure title
        plt.rc('axes.formatter', limits=(-3, 3))

        # Z parf poss & fits
        ie_budget_plot_z(z_parf_poss_samp + z_parf_fits_samp, 'b.', 0.75, row=5, col=6, sub=1,
                         ylab='Z parf poss & fits ($\mu$m)',
                         ylim=[-1.5e6 * (z_parf_poss + z_parf_fits), 1.5e6 * (z_parf_poss + z_parf_fits)])
        # Z aiv
        ie_budget_plot_z(z_aiv_samp, 'b.', 0.75, row=5, col=6, sub=2,
                         ylab='Z AIV ($\mu$m)', ylim=[-1.5e6 * z_aiv, 1.5e6 * z_aiv])
         # Z lookup
        ie_budget_plot_z(z_lookup_samp, 'b.', 0.75, row=5, col=6, sub=3,
                         ylab='Z lookup ($\mu$m)', ylim=[-1.5e6 * z_lookup, 1.5e6 * z_lookup])
        # Z defoc poss
        ie_budget_plot_z(z_defoc_poss_samp, 'b.', 0.75, row=5, col=6, sub=4,
                         ylab='Z PosS ($\mu$m)', ylim=[-50, 150])
        # Z grav
        ie_budget_plot_z(z_grav_samp[:,0], 'g.', 0.5, row=5, col=6, sub=7,
                         ylab='Z grav ($\mu$m)', ylim=[-1.5e6 * z_grav, 1.5e6 * z_grav])
        # Z pfhs
        ie_budget_plot_z(z_pfhs_samp[:,0], 'g.', 0.5, row=5, col=6, sub=8,
                         ylab='Z PFHS ($\mu$m)', ylim=[-1.5e6 * z_pfhs, 1.5e6 * z_pfhs])
        # Z therm
        ie_budget_plot_z(z_therm_samp[:,0], 'g.', 0.5, row=5, col=6, sub=9,
                         ylab='Z therm ($\mu$m)', ylim=[-1.5e6 * z_therm, 1.5e6 * z_therm])
        # Z InRo
        ie_budget_plot_z(z_inro_samp[:,0], 'g.', 0.5, row=5, col=6, sub=10,
                         ylab='Z InRo ($\mu$m)', ylim=[-1.5e6 * z_inro, 1.5e6 * z_inro])

        # XY tgt
        ie_budget_plot_xy(x_tgtcoo_samp, y_tgtcoo_samp, 'b.', 0.75, row=5, col=6, sub=13,
                          xlab='X tgt coo ($\mu$m)', ylab='Y tgt coo ($\mu$m)',
                          xlim=[-1.5e6 * r_tgtcoo, 1.5e6 * r_tgtcoo], ylim=[-1.5e6 * r_tgtcoo, 1.5e6 * r_tgtcoo])
        # XY skymap
        ie_budget_plot_xy(x_skymap_samp, y_skymap_samp, 'b.', 0.75, row=5, col=6, sub=14,
                          xlab='X skymap ($\mu$m)', ylab='Y skymap ($\mu$m)',
                          xlim=[-1.5e6 * r_skymap, 1.5e6 * r_skymap], ylim=[-1.5e6 * r_skymap, 1.5e6 * r_skymap])
        # XY PosS
        ie_budget_plot_xy(x_poss_samp, y_poss_samp, 'b.', 0.75, row=5, col=6, sub=15,
                          xlab='X PosS ($\mu$m)', ylab='Y PosS ($\mu$m)',
                          xlim=[-1.5e6 * r_poss_aao[1], 1.5e6 * r_poss_aao[1]], ylim=[-1.5e6 * r_poss_aao[1], 1.5e6 * r_poss_aao[1]])
        # XY therm
        ie_budget_plot_xy(x_therm_samp[:,0], y_therm_samp[:,0], 'g.', 0.5, row=5, col=6, sub=19,
                          xlab='X therm ($\mu$m)', ylab='Y therm ($\mu$m)',
                          xlim=[-1.5e6 * r_therm, 1.5e6 * r_therm], ylim=[-1.5e6 * r_therm, 1.5e6 * r_therm])
        # XY vibr
        ie_budget_plot_xy(x_vibr_samp[:,0], y_vibr_samp[:,0], 'g.', 0.5, row=5, col=6, sub=20,
                          xlab='X vibr ($\mu$m)', ylab='Y vibr ($\mu$m)',
                          xlim=[-1.5e6 * r_vibr, 1.5e6 * r_vibr], ylim=[-1.5e6 * r_vibr, 1.5e6 * r_vibr])
        # XY guid
        ie_budget_plot_xy(x_guid_samp[:,0], y_guid_samp[:,0], 'g.', 0.5, row=5, col=6, sub=21,
                          xlab='X guid ($\mu$m)', ylab='Y guid ($\mu$m)',
                          xlim=[-2e6 * r_guid, 2e6 * r_guid], ylim=[-2e6 * r_guid, 2e6 * r_guid])
        # XY InRo
        ie_budget_plot_xy(x_inro_samp[:,0], y_inro_samp[:,0], 'g.', 0.5, row=5, col=6, sub=22,
                          xlab='X InRo ($\mu$m)', ylab='Y InRo ($\mu$m)',
                          xlim=[-2e6 * r_inro, 2e6 * r_inro], ylim=[-2e6 * r_inro, 2e6 * r_inro])
        # XY ADR
        ie_budget_plot_xy(x_adr_samp[:,0], y_adr_samp[:,0], 'g.', 0.5, row=5, col=6, sub=23,
                          xlab='X ADR ($\mu$m)', ylab='Y ADR ($\mu$m)',
                          xlim=[-1.5e6 * r_adr, 1.5e6 * r_adr], ylim=[-1.5e6 * r_adr, 1.5e6 * r_adr])
        # XY plate
        ie_budget_plot_xy(x_plate_samp[:,0], y_plate_samp[:,0], 'g.', 0.5, row=5, col=6, sub=24,
                          xlab='X plate ($\mu$m)', ylab='Y plate ($\mu$m)',
                          xlim=[-1.5e6 * r_plate, 1.5e6 * r_plate], ylim=[-1.5e6 * r_plate, 1.5e6 * r_plate])

        # Total offsets
        #  X, Y
        ie_budget_plot_xy(x_off, y_off, 'k.', 0.01, row=5, col=6, sub=25,
                          xlab='X offset ($\mu$m)', ylab='Y offset ($\mu$m)',
                          xlim=[-50, 50], ylim=[-50, 50])
        # X, Z
        ie_budget_plot_xy(x_off, z_off, 'k.', 0.01, row=5, col=6, sub=26,
                          xlab='X offset ($\mu$m)', ylab='Z offset ($\mu$m)',
                          xlim=[-50, 50], ylim=[-150, 200])
        # Y, Z
        ie_budget_plot_xy(y_off, z_off, 'k.', 0.01, row=5, col=6, sub=27,
                          xlab='Y offset ($\mu$m)', ylab='Z offset ($\mu$m)',
                          xlim=[-50, 50], ylim=[-150, 200])
        # adjust space in between plots
        plt.subplots_adjust(top=0.99, bottom=0.08, left=0.08, right=0.99, hspace=0.35,
                            wspace=0.45)
        # save figure
        plt.savefig('results/simu/' + self.design + '_xyz_offsets_' + self.spectro + '_fib' + str(self.fibdiam) +
                    '_Z' + self.zenith + '_Field_' + self.field + '_Seeing_' + str(self.seeing) + '.png')
        plt.close()
#        plt.rcdefaults()

        # --- Plot Histograms ---

        plt.figure(figsize=(10, 8))
        plt.subplot(5, 6, 1)
        plt.hist(np.ravel(z_parf_poss_samp * 1e6), 30, color='b', histtype='step')
        plt.xlabel('Z parf poss ($\mu$m)')
        plt.subplot(5, 6, 2)
        plt.hist(np.ravel(z_parf_fits_samp * 1e6), 30, color='b', histtype='step')
        plt.xlabel('Z parf fits ($\mu$m)')
        plt.subplot(5, 6, 3)
        plt.hist(np.ravel(z_lookup_samp * 1e6), 30, color='b', histtype='step')
        plt.xlabel('Z lookup ($\mu$m)')
        plt.subplot(5, 6, 4)
        plt.hist(np.ravel(z_defoc_poss_samp * 1e6), 30, color='b', histtype='step')
        plt.xlabel('Z PosS ($\mu$m)')
        plt.subplot(5, 6, 7)
        plt.hist(np.ravel(z_grav_samp * 1e6), 30, color='g', histtype='step')
        plt.xlabel('Z grav ($\mu$m)')
        plt.subplot(5, 6, 8)
        plt.hist(np.ravel(z_pfhs_samp * 1e6), 30, color='g', histtype='step')
        plt.xlabel('Z PFHS ($\mu$m)')
        plt.subplot(5, 6, 9)
        plt.hist(np.ravel(z_therm_samp * 1e6), 30, color='g', histtype='step')
        plt.xlabel('Z therm ($\mu$m)')
        plt.subplot(5, 6, 10)
        plt.hist(np.ravel(z_inro_samp * 1e6), 30, color='g', histtype='step')
        plt.xlabel('Z InRo ($\mu$m)')

        plt.subplot(5, 6, 13)
        plt.hist(np.ravel(x_tgtcoo_samp * 1e6), 30, color='b', histtype='step')
        plt.hist(np.ravel(y_tgtcoo_samp * 1e6), 30, color='b', alpha=0.5, histtype='step')
        plt.xlabel('XY tgt coo ($\mu$m)')
        plt.subplot(5, 6, 14)
        plt.hist(np.ravel(x_skymap_samp * 1e6), 30, color='b', histtype='step')
        plt.hist(np.ravel(y_skymap_samp * 1e6), 30, color='b', alpha=0.5, histtype='step')
        plt.xlabel('XY skymap ($\mu$m)')
        plt.subplot(5, 6, 15)
        plt.hist(np.ravel(x_poss_samp * 1e6), 30, color='b', histtype='step')
        plt.hist(np.ravel(y_poss_samp * 1e6), 30, color='b', alpha=0.5, histtype='step')
        plt.xlabel('XY PosS ($\mu$m)')

        plt.subplot(5, 6, 19)
        plt.hist(np.ravel(x_therm_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_therm_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY therm ($\mu$m)')
        plt.subplot(5, 6, 20)
        plt.hist(np.ravel(x_vibr_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_vibr_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY vibr ($\mu$m)')
        plt.subplot(5, 6, 21)
        plt.hist(np.ravel(x_poss_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_poss_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY guid ($\mu$m)')
        plt.subplot(5, 6, 22)
        plt.hist(np.ravel(x_inro_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_inro_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY InRo ($\mu$m)')
        plt.subplot(5, 6, 23)
        plt.hist(np.ravel(x_adr_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_adr_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY ADR ($\mu$m)')
        plt.subplot(5, 6, 24)
        plt.hist(np.ravel(x_plate_samp * 1e6), 30, color='g', histtype='step')
        plt.hist(np.ravel(y_plate_samp * 1e6), 30, color='g', alpha=0.5, histtype='step')
        plt.xlabel('XY plate ($\mu$m)')

        plt.subplot(5, 6, 25)
        plt.hist(np.ravel(x_off * 1e6), 30, color='k', histtype='step')
        plt.xlabel('X offset ($\mu$m)')
        plt.subplot(5, 6, 26)
        plt.hist(np.ravel(y_off * 1e6), 30, color='k', histtype='step')
        plt.xlabel('Y offset ($\mu$m)')
        plt.subplot(5, 6, 27)
        plt.hist(np.ravel(z_off * 1e6), 30, color='k', histtype='step')
        plt.xlabel('Z offset ($\mu$m)')
        # adjust space in between plots
        plt.subplots_adjust(top=0.99, bottom=0.08, left=0.08, right=0.99, hspace=0.35,
                            wspace=0.45)
        plt.savefig('results/simu/' + self.design + '_xyz_offsets_hist_' + self.spectro + '_fib' + str(self.fibdiam) +
                    '_Z' + self.zenith + '_Field_' + self.field + '_Seeing_' + str(self.seeing) + '.png')

        plt.close()

        # save offsets
        self.x_off = x_off
        self.y_off = y_off
        self.z_off = z_off

    def compute_ie(self, optim_ie):
        """ Compute the Injection Efficiency given the random values picked up
        """

        # Open FITS file
        hdu = fits.open(self.cubeconvfile)
        cubeconv = hdu[0].data
        hdr = hdu[0].header
        # Get optimal position
        optim_pos = [hdr['OPTX' + self.spectro[0] + str(self.fibdiam)[0] + str(self.fibdiam)[2:4]],
                     hdr['OPTY' + self.spectro[0] + str(self.fibdiam)[0] + str(self.fibdiam)[2:4]]]

        # Create all the apertures for that simulation
        t0 = time.time()

        aper_list = [phot.CircularAperture([(optim_pos[0] + self.x_off[i, j]/pix_scale,
                                             optim_pos[1] + self.y_off[i, j]/pix_scale)
                                            for i in range(nstep)], r=self.fibdiam / 2 * plate_scale / pix_scale)
                     for j in range(nsimu)]

        t1 = time.time()

        # Perform aperture photometry at all defocus values (and all wavelengths) -- normalizing properly!
        ie_all = [[[phot.aperture_photometry(cubeconv[f, w, :, :] / np.sum(self.cube[f][w, :, :]), aper_list[j])
                    for j in range(nsimu)]
                   for f in range(nfoci)]
                  for w in range(nwave)]
        # This gives 5x11 "tables" of aperture photometry object, each table being nsimu*nstep values

        t2 = time.time()

        # Interpol at correct defocus
        ie_def = [[[np.interp(self.z_off[i, j], list(foci.values()),
                              [ie_all[w][0][j]['aperture_sum'].quantity[i],
                               ie_all[w][1][j]['aperture_sum'].quantity[i],
                               ie_all[w][2][j]['aperture_sum'].quantity[i],
                               ie_all[w][3][j]['aperture_sum'].quantity[i],
                               ie_all[w][4][j]['aperture_sum'].quantity[i]])
                    for i in range(nstep)]
                   for j in range(nsimu)]
                  for w in range(nwave)]

        t3 = time.time()

        # Average all steps
        ie_step = [[np.mean(ie_def[w][j])
                    for j in range(nsimu)]
                   for w in range(nwave)]

        t4 = time.time()

        # Average all simu
        ie_simu = [np.mean(ie_step[w])
                   for w in range(nwave)]

        ie_simu_stddev = [np.std(ie_step[w])
                          for w in range(nwave)]

        t5 = time.time()

        # Plot IE
        fig = plt.figure(figsize=(10, 6))
        # IE curves
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.plot(wave, ie_simu, 'r', zorder=10, linewidth=2, label='Overall mean')
        ax1.plot(wave, optim_ie, 'g-', zorder=10, linewidth=2, label='Optimal position')
        for i in range(nsimu):
            ax1.plot(wave, [tmp[i] for tmp in ie_step], alpha=.25, color='b', zorder=5,
                     label='All simulations' if i == 0 else '')
            for j in range(nstep):
                ax1.plot(wave, [tmp[i][j] for tmp in ie_def], alpha=.1, color='k', zorder=1,
                         label='All samples' if (i, j) == (0, 0) else '')
        ax1.set_xlabel('Wavelengths (nm)')
        ax1.set_ylabel('Injection efficiency')
        ax1.legend()
        # Standard deviation (absolute)
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(wave, [np.std(ie_step[i]) for i in range(11)], 'r', zorder=10, label='Between all simulations')
        for i in range(nsimu):
            ax2.plot(wave, [np.std(ie_def[j][i]) for j in range(11)], alpha=.25, color='b', zorder=5,
                     label='Between all samples' if i == 0 else '')
        ax2.set_xlabel('Wavelengths (nm)')
        ax2.set_ylabel('StdDev (injection efficiency)')
        ax2.legend()
        # Save plot
        fig.savefig('results/simu/' + self.design + '_ie_simu_' + self.spectro + '_fib' + str(self.fibdiam) +
                    '_Z' + self.zenith + '_Field_' + self.field + '_Seeing_' + str(self.seeing) + '.png')
        plt.close(fig)

        t6 = time.time()

        print(" -- breakdown of computing IE:", t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, t6 - t5)

        # Save IE
        self.ie_simu = ie_simu
        self.ie_simu_stddev = ie_simu_stddev

    def simu_offset_adr(self):

        # distance from center in meters
        field_dist = float(self.field[2:]) * plate_scale * 3600.

        # also need axis to get orientation
        field_axis = self.field[0]
        theta_ref = np.pi / 2.
        if field_axis == 'X':
            theta_ref = 0.

        # random variation at edge of field
        edge = np.random.random_sample(nstep) * r_adr * 2 - r_adr

        # sort them to account for smoothness of variations
        edge = np.sort(edge)

        # scale to FIELD_DIST
        roff = field_dist * edge / fov_rad_m

        # centripetal or centrifugal?
        if field_axis == 'Y':
            roff *= -1

        # get X and Y offset
        xoff = roff * np.cos(theta_ref)
        yoff = roff * np.sin(theta_ref)

        return xoff, yoff

    def simu_offset_plate(self):

        # Assuming constant plate-scale to make it simpler

        # distance from center in meters
        field_dist = float(self.field[2:]) * plate_scale * 3600.

        # also need axis to get orientation
        field_axis = self.field[0]
        theta_ref = np.pi / 2.
        if field_axis == 'X':
            theta_ref = 0.

        # random variation at edge of field
        edge = np.random.random_sample(nstep) * r_plate * 2 - r_plate

        # sort them to account for smoothness of variations
        edge = np.sort(edge)
        flip = np.round(np.random.random_sample(1))  # flip a coin to increase/decrease
        if flip:
            edge = edge[::-1]

        # scale to FIELD_DIST
        roff = field_dist * edge / fov_rad_m

        # get X and Y offset
        xoff = roff * np.cos(theta_ref)
        yoff = roff * np.sin(theta_ref)

        return xoff, yoff

    def simu_offset_inro(self):
        """ Computes offset due to InRo
        """
        # distance from center in meters
        field_dist = float(self.field[2:]) * plate_scale * 3600

        # also need axis to get orientation of InRo effect
        field_axis = self.field[0]
        theta_ref = np.pi / 2.
        if field_axis == 'X':
            theta_ref = 0.

        # radii values (constant) in microns
        rr = np.ones(nstep) * field_dist

        # theta values, near theta_ref, in radians
        if field_dist == 0:
            theta = rr
        else:
            theta = theta_ref + np.random.normal(size=nstep) * np.arctan(r_inro/field_dist)

        # convert to X and Y, in microns
        xoff = rr * (np.cos(theta) - np.cos(theta_ref))
        yoff = rr * (np.sin(theta) - np.sin(theta_ref))

        return xoff, yoff

    def simu_defoc_poss(self):

        # Fiber pitch distribution
        #   get all results from allocation efficiency simulations
        path = '/Users/nflagey/WORK/PROGRAM/IDL/eklib/mse/'
        files = os.listdir(path)
        all_d = []
        if self.spectro == 'MR' or self.spectro == 'LR':
            root = 'LR'
        else:
            root = 'HR'
        for name in files:
            if fnmatch.fnmatch(name, 'mse_alloc_simu_aao*'+root+'*save'):
                result = os.path.join(path, name)
                alloc = scipy.io.readsav(result)
                all_d += list(alloc['all_d'])
        all_d = np.asarray(all_d) * plate_scale  # arcsec --> m

        # Randomly pick pitch values
        pitch = np.random.choice(all_d, size=nsimu)
        # Compute defocus
        depth = spine_length - np.sqrt(spine_length**2 - pitch**2)

        # Find median value to define focus=0
        pitch_med = np.median(all_d)
        depth_med = spine_length - np.sqrt(spine_length**2 - pitch_med**2)
  
        # Remove defocus at median pitch
        defoc = depth - depth_med

        return defoc

    def simu_defoc_inrotilt(self):

        # distance from center in meters
        field_dist = float(self.field[2:]) * plate_scale * 3600

        # pick a random starting angle
        theta_init = np.random.random_sample(1) * 360.  # in degrees

        # angle range
        # - restore mse_keyhole results
        keyhole = scipy.io.readsav('/Users/nflagey/WORK/MSE/MSE_Keyhole/mse_keyhole_for_injection_efficiency.sav')
        # - get InRo range for ALL_ZD_MEAN within 5 degrees of input ZD (60mn)
        inro_60 = keyhole['all_rot'][1, 1, :, :, :]  # starting point above ZD60 because we do not observe below
        zd_mean_60 = keyhole['all_zd_mean'][1, 1, :, :, :]
        inro_range_60 = inro_60[(keyhole['all_alt1'] >= 0.) & (abs(zd_mean_60 - float(self.zenith)) < 5.)]

        # pick random InRo ranges within that distribution
        theta_range_60 = np.random.choice(inro_range_60, size=1)

        # for each simulation (1 per observation)
        # now we need to define the sample size for the observation
        theta_step_60 = theta_range_60 / (nstep - 1)
        # now we define all THETA values
        theta_60 = theta_init + theta_step_60 * np.arange(nstep)
        # now that we know all the THETA values, we can compute the defocus values
        defoc_60 = z_inro * field_dist / fov_rad_m * np.cos(theta_60 / 180 * np.pi)

        return defoc_60


def ie_budget_plot_xy(x, y, mark, alpha=0.25, ls='None', row=4, col=5, sub=1,
                      xlab='', ylab='', xlim=[0., 1.], ylim=[0., 1.]):
    plt.subplot(row, col, sub)
    plt.plot(x * 1e6, y * 1e6, mark, alpha=alpha, ls=ls, markersize=2)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(xlim)
    plt.ylim(ylim)


def ie_budget_plot_z(z, mark, alpha=0.25, ls='None', row=4, col=5, sub=1, ylab='', ylim=[0., 1.]):
    plt.subplot(row, col, sub)
    plt.plot(z * 1e6, mark, alpha=alpha, ls=ls, markersize=2)
    plt.ylabel(ylab)
    plt.xlabel('Simulations')
    plt.ylim(ylim)

def ie_moffat(seeing_wave):
    # Moffat parameters
    gamma = 3.
    hwhm = seeing_wave * plate_scale / pix_scale / 2.
    alpha = hwhm / np.sqrt(2. ** (1. / gamma) - 1)
    size = int(8. * np.round(alpha) + 1.)
    # Generate kernel from scratch (astropy models leads to very peaky kernel)
    kernel = [[(gamma - 1) / (np.pi * alpha ** 2) *
               (1 + ((i - (size - 1) / 2) ** 2 + (j - (size - 1) / 2) ** 2) / alpha ** 2) ** (- gamma)
               for i in range(size)] for j in range(size)]
    return kernel
