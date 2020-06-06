#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy import constants as const
from astropy.io import fits
from astropy.table import Table
from scipy import signal
import copy
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot
from bokeh.embed import components
from bokeh.models.ranges import Range1d
from bokeh.models.axes import LogAxis


class MseSpectrum:

    def __init__(self, sessionID=-1, tgtmag=24., band='g', template='flat', src_type='point', redshift=0., airmass=1.0,
                 skymag=20.7, seeing=0.4, coating='ZeCoat', fibdiam=1.0, spectro='LR', specbin=1., spatbin=2.,
                 meth='getSNR', snr=2., etime=3600., badtgtmag=0.):
        # SessionID and path (for local uses sessionID=-1, otherwise sessionID is set by the HTML)
        self.sessionID = sessionID
        global path
        path = './'
        if self.sessionID == -1:
            path = '/mse_itc_dev/'
        # Grids
        self.wgrid = np.array([])
        self.reselgrid = np.array([])
        self.egrid = np.array([])
        self.armgrid = np.array([])
        # Target
        self.tgtmag = tgtmag
        self.band = band
        self.template = template
        self.src_type = src_type
        self.redshift = redshift
        self.tgtflux = np.array([])
        self.tgtcount = np.array([])
        self.tgtdetec = np.array([])
        self.tgtnoise = np.array([])
        # Sky
        self.airmass = airmass
        self.skymag = skymag
        self.seeing = seeing
        self.skyflux = np.array([])
        self.skycount = np.array([])
        self.skytrans = np.array([])
        self.skydetec = np.array([])
        self.skynoise = np.array([])
        # Telescope
        self.coating = coating
        self.fibdiam = fibdiam
        self.spectro = spectro
        self.dark = np.array([])
        self.darknoise = np.array([])
        self.thermal = np.array([])
        self.thermalnoise = np.array([])
        self.readout = np.array([])
        self.specbin = specbin
        self.spatbin = spatbin
        # Throughput and other curves
        self.thr_struc = np.array([])
        self.thr_m1 = np.array([])
        self.thr_pfhs = np.array([])
        self.thr_pfue = np.array([])
        self.thr_poss = np.array([])
        self.thr_fiber = np.array([])
        self.thr_spectro = np.array([])
        self.inj = np.array([])
        # Noise that depends on other sources (0 in default use of online ITC)
        self.badtgtmag = badtgtmag
        self.xtalk = np.array([])
        self.ghost = np.array([])
        self.teldiffuse = np.array([])
        self.instdiffuse = np.array([])

        # Methode
        self.meth = meth
        self.snr = snr
        self.etime = etime

        # Create grids
        self.create_grids()
        # Create target spectrum
        self.create_target()
        # Create skyspectrum
        self.create_sky()


    def create_grids(self):
        """Makes the wavelength grid for each spectrograph.
        """

        # Read CSV files (digitize curves from CoDR report)
        # Depends on spectrograph
        if self.spectro == 'LR':
            # wav in micron
            blue = Table.read(path + 'THROUGHPUT/LMR_resol_LRblue.csv', format='csv')
            green = Table.read(path + 'THROUGHPUT/LMR_resol_LRgreen.csv', format='csv')
            red = Table.read(path + 'THROUGHPUT/LMR_resol_LRred.csv', format='csv')
            nirJ = Table.read(path + 'THROUGHPUT/LMR_resol_LRj.csv', format='csv')
            nirH = Table.read(path + 'THROUGHPUT/LMR_resol_MRh.csv', format='csv')
            # wav in angstroms
            blue['wav'] *= 1e4
            green['wav'] *= 1e4
            red['wav'] *= 1e4
            nirJ['wav'] *= 1e4
            nirH['wav'] *= 1e4
        elif self.spectro == 'MR':
            # wav in micron
            blue = Table.read(path + 'THROUGHPUT/LMR_resol_MRblue.csv', format='csv')
            green = Table.read(path + 'THROUGHPUT/LMR_resol_MRgreen.csv', format='csv')
            red = Table.read(path + 'THROUGHPUT/LMR_resol_MRred.csv', format='csv')
            # wav in angstroms
            blue['wav'] *= 1e4
            green['wav'] *= 1e4
            red['wav'] *= 1e4
        else:
            # wav in nm
            blue = Table.read(path + 'THROUGHPUT/HR_resol_blue.csv', format='csv')
            green = Table.read(path + 'THROUGHPUT/HR_resol_green.csv', format='csv')
            red = Table.read(path + 'THROUGHPUT/HR_resol_red.csv', format='csv')
            # wav in angstroms
            blue['wav'] *= 10.
            green['wav'] *= 10.
            red['wav'] *= 10.

        # Compute resolution element grid, in angstrom
        reselgrid_blue = blue['wav'] / blue['resol']
        reselgrid_green = green['wav'] / green['resol']
        reselgrid_red = red['wav'] / red['resol']
        if self.spectro == 'LR':
            reselgrid_nir_j = nirJ['wav'] / nirJ['resol']
            reselgrid_nir_h = nirH['wav'] / nirH['resol']

        # Sampling grid is 1/2 of resolution element grid
        sampgrid_blue = reselgrid_blue / 2.
        sampgrid_green = reselgrid_green / 2.
        sampgrid_red = reselgrid_red / 2.
        if self.spectro == 'LR':
            sampgrid_nir_j = reselgrid_nir_j / 2.
            sampgrid_nir_h = reselgrid_nir_h / 2.

        # Now we can built the wavelength grid
        # blue
        wgrid_blue = np.array([np.min(blue['wav'])])
        while True:
            wgrid_blue = np.append(wgrid_blue, wgrid_blue[-1] + np.interp(wgrid_blue[-1], blue['wav'], sampgrid_blue))
            if wgrid_blue[-1] > np.max(blue['wav']):
                break
        # green
        wgrid_green = np.array([np.min(green['wav'])])
        while True:
            wgrid_green = np.append(wgrid_green, wgrid_green[-1] + np.interp(wgrid_green[-1], green['wav'], sampgrid_green))
            if wgrid_green[-1] > np.max(green['wav']):
                break
        # red
        wgrid_red = np.array([np.min(red['wav'])])
        while True:
            wgrid_red = np.append(wgrid_red, wgrid_red[-1] + np.interp(wgrid_red[-1], red['wav'], sampgrid_red))
            if wgrid_red[-1] > np.max(red['wav']):
                break
        # NIR
        if self.spectro == 'LR':
            wgrid_nir_j = np.array([np.min(nirJ['wav'])])
            while True:
                wgrid_nir_j = np.append(wgrid_nir_j, wgrid_nir_j[-1] + np.interp(wgrid_nir_j[-1], nirJ['wav'], sampgrid_nir_j))
                if wgrid_nir_j[-1] > np.max(nirJ['wav']):
                    break
            wgrid_nir_h = np.array([np.min(nirH['wav'])])
            while True:
                wgrid_nir_h = np.append(wgrid_nir_h, wgrid_nir_h[-1] + np.interp(wgrid_nir_h[-1], nirH['wav'], sampgrid_nir_h))
                if wgrid_nir_h[-1] > np.max(nirH['wav']):
                    break

        #  Interpolate resolution element grid and convert to angstroms
        reselgrid_blue = np.interp(wgrid_blue, blue['wav'], reselgrid_blue)
        reselgrid_green = np.interp(wgrid_green, green['wav'], reselgrid_green)
        reselgrid_red = np.interp(wgrid_red, red['wav'], reselgrid_red)
        if self.spectro == 'LR':
            reselgrid_nir_j = np.interp(wgrid_nir_j, nirJ['wav'], reselgrid_nir_j)
            reselgrid_nir_h = np.interp(wgrid_nir_h, nirH['wav'], reselgrid_nir_h)
        
        # Combine arms (convert wavelength into angstrom)
        self.wgrid = np.append(np.append(wgrid_blue, wgrid_green), wgrid_red)
        self.reselgrid = np.append(np.append(reselgrid_blue, reselgrid_green), reselgrid_red)
        self.armgrid = np.append(np.append(wgrid_blue * 0, wgrid_green * 0 + 1), wgrid_red * 0 + 2)
        if self.spectro == 'LR':
            self.wgrid = np.append(np.append(self.wgrid, wgrid_nir_j), wgrid_nir_h)
            self.reselgrid = np.append(np.append(self.reselgrid, reselgrid_nir_j), reselgrid_nir_h)
            # Create a variable that knows the arm
            self.armgrid = np.append(np.append(self.armgrid, wgrid_nir_j*0+3), wgrid_nir_h*0+4)

        # Sort wavelength and apply to all grids
        self.egrid = 1e17 * const.h.value * const.c.value / self.wgrid
        s = np.argsort(self.wgrid)
        self.wgrid = self.wgrid[s]
        self.reselgrid = self.reselgrid[s]
        self.egrid = self.egrid[s]
        self.armgrid = self.armgrid[s]


    def create_target(self):
        """Create target spectrum.
        """

        # Open filter and interpolate onto nominal wavelength grid
        filt = Table.read(path + 'FILTERS/' + self.band + '.dat', format='ascii')
        filtrans = np.interp(self.wgrid, filt['col1'], filt['col2'])

        # Pick up correct template
        if self.template == 'flat':
            tgttemp_lam = self.wgrid.copy()
            tgttemp_flux, tgttemp_count = mag2flux(tgttemp_lam, self.tgtmag)  # in erg/s/cm2/A and ph/s/cm2/A
        else:
            temp = Table.read(path + 'TEMPLATE2/' + self.template + '.dat', format='ascii')
            if self.template == 'qso1':
                tgttemp_lam = temp['col1'] * 1e4  # microns to Angstroms
                tgttemp_flux = temp['col2'] * 1e10 * const.c.value / tgttemp_lam ** 2.  # f_nu to f_lambda
            elif self.template == 'qso2' or self.template == 'elliptical' or self.template == 'spiral_sc' or self.template == 'HII' or self.template == 'PN':
                tgttemp_lam = temp['col1'] * 10.  # nm to Angstroms
                tgttemp_flux = temp['col2']
            else:
                tgttemp_lam = temp['col1']
                tgttemp_flux = temp['col2']

        # Apply redshift (assuming flux is given as observed, so no conversion)
        tgttemp_lam *= 1 + self.redshift

        # Interpolate onto nominal wavelength grid (put to 0 if target spectrum is narrower)
        tgttemp_flux = np.interp(self.wgrid, tgttemp_lam, tgttemp_flux, left=0, right=0)
        tgttemp_count = tgttemp_flux / self.egrid

        # Find total counts for input magnitude
        tgtmag_flux, tgtmag_count = mag2flux(self.wgrid, self.tgtmag)
        tgtmag_count_sum = np.sum(tgtmag_count * filtrans)
        # Find total counts for template
        tgttemp_count_sum = np.sum(tgttemp_count * filtrans)

        # Scale template to match magnitude
        self.tgtflux = tgttemp_flux / tgttemp_count_sum * tgtmag_count_sum
        self.tgtcount = self.tgtflux * self.egrid

    def create_sky(self):
        """Create sky spectrum.
        """
        # Open V filter
        vfilt = Table.read(path + 'FILTERS/Vbm.raw', format='ascii')

        # Open correct sky file
        am = str(int(self.airmass * 10.))
        sky = str(int(self.skymag * 10.))
        hdul = fits.open(path + 'SKY/skytable_am' + am + '_sky' + sky + '.fits')
        data = hdul[1].data
        skytemp_lam = data['lam'] * 1e4  # in microns --> Angstroms
        skytemp_count = data['flux'] / 1e8  # in ph/s/m2/um/arcsec2 --> ph/s/cm2/A/arcsec2
        skytemp_trans = data['trans']  # in 0 to 1

        # Need to work on template wavgrid to avoid interpolation issues with the template during normalization
        # Interpolate filter onto sky grid for normalization
        vfiltrans = np.interp(skytemp_lam, vfilt['col1'], vfilt['col2'])
        # Find total counts for template
        skytemp_count_sum = np.trapz(skytemp_count * vfiltrans, skytemp_lam)
        # Find total flux (ergs/s/cm2/A) and counts for input magnitude
        skymag_flux, skymag_count = mag2flux(skytemp_lam, self.skymag)
        skymag_count_sum = np.trapz(skymag_count * vfiltrans, skytemp_lam)

        # Now we can interpolate onto correct wavelength grid
        skytemp_trans = np.interp(self.wgrid, skytemp_lam, skytemp_trans)
        skytemp_count = np.interp(self.wgrid, skytemp_lam, skytemp_count) / skytemp_count_sum * skymag_count_sum

        # Scale template to match magnitude
        self.skycount = skytemp_count
        self.skyflux = self.skycount * self.egrid  # ergs/s/cm2/A
        self.skytrans = skytemp_trans

    def apply_atmos_ext(self):
        """Apply atmospheric extinction.
        """
        # Make a copy of the input spectrum so it's not modified
        spec = copy.deepcopy(self)

        # Apply extinction
        spec.tgtflux *= spec.skytrans

        return spec

    def apply_throughput_front(self):
        """Apply throughput of TEL, M1, and PFUE.
        """

        # Make a copy of the input spectrum so it's not modified
        spec = copy.deepcopy(self)

        # Open throughput files for TEL, M1, and PFUE
        struc = Table.read(path + 'THROUGHPUT/mse_itc_throughput_struc.dat', format='ascii')
        if self.coating == 'ZeCoat':
            m1 = Table.read(path + 'THROUGHPUT/mse_itc_throughput_m1_zecoat.dat', format='ascii')
        else:
            m1 = Table.read(path + 'THROUGHPUT/mse_itc_throughput_m1_gemini.dat', format='ascii')
        pfue = Table.read(path + 'THROUGHPUT/mse_itc_throughput_pfue_WS.dat', format='ascii')
        pfhs = Table.read(path + 'THROUGHPUT/mse_itc_throughput_pfhs.dat', format='ascii')
        # Interpolate
        thr_struc = np.interp(spec.wgrid, struc['lamA'], struc['thr'])
        thr_m1 = np.interp(spec.wgrid, m1['lamA'], m1['thr'])
        thr_pfhs = np.interp(spec.wgrid, pfhs['lamA'], pfhs['thr'])
        thr_pfue = np.interp(spec.wgrid, pfue['Wavelength'] * 10000., pfue['Transmission'] * pfue['Reflexion'] * pfue['Vignetting'])

        # Apply throughput to sky and target
        spec.tgtflux *= thr_struc * thr_m1 * thr_pfhs * thr_pfue
        spec.skyflux *= thr_struc * thr_m1 * thr_pfhs * thr_pfue

        # Save for debug
        spec.thr_struc = thr_struc
        spec.thr_m1 = thr_m1
        spec.thr_pfue = thr_pfue
        spec.thr_pfhs = thr_pfhs

        return spec

    def apply_injeff(self):
        """Applies the injection efficiency.
        """

        # Make a copy of the input spectrum so it's not modified
        spec = copy.deepcopy(self)

        # Load IE curves from file (dictionary of dictionaries of dictionaries ...)
        ie_wav = np.array([360., 370., 400., 445., 551., 658., 806., 1000., 1214., 1477., 1784.])
        ie_for_itc = np.load(path + 'THROUGHPUT/no_segments_injeff_curve_for_itc_all.npy').flat[0]
        # ie_avg_all['iq0.45']['ZD30']['LR1.00']

        # Selecting correct curves for this observation
        # Fiber diameter & spectro
        specfib_str = spec.spectro + "{:.2f}".format(spec.fibdiam)
        # ZD/Airmass
        zd_airmass = {1.0: '00', 1.2: '30', 1.5: '50', 2.0: '60'}
        zd_str = "ZD" + zd_airmass[spec.airmass]
        # get all IQs in dictionary
        iqs = sorted(list(ie_for_itc.keys()))
        iqvals = [float(iqs[i][2:]) for i in range(len(iqs))]
        # isolate IE curves of interest
        ie_itc = {}
        for iq in iqs:
            ie_itc[iq] = ie_for_itc[iq][zd_str][specfib_str]

        # Interpolate at desired IQ
        ie_coarse = [np.interp(spec.seeing, iqvals, [ie_itc[iq][w] for iq in iqs]) for w in range(len(ie_wav))]

        # Interpolate over correct wavelength grid
        ie_curve = np.interp(spec.wgrid, ie_wav * 10., ie_coarse)

        # Apply injection efficiency to target
        spec.tgtflux *= ie_curve
        # Apply injection efficiency to sky (only apply fiber area)
        spec.skyflux *= np.pi * (spec.fibdiam / 2.) ** 2.

        # Save for verification
        spec.inj = ie_curve

        return spec


    def apply_throughput_back(self):
        """Apply throughput of PosS (no Inj.Eff.), FiTS, and Spectro.
        """

        # Make a copy of the input spectrum so it's not modified
        spec = copy.deepcopy(self)

        # Open throughput files for FiTS, Spectro, and Poss
        poss = Table.read(path + 'THROUGHPUT/mse_itc_throughput_poss.dat', format='ascii')
        if spec.spectro == 'LR' or spec.spectro == 'MR':
            fiber = Table.read(path + 'THROUGHPUT/mse_itc_throughput_fits_lmr.dat', format='ascii')
            spectro = Table.read(path + 'THROUGHPUT/mse_itc_throughput_spec_lr.dat', format='ascii')
            if spec.spectro == 'MR':
                spectro = Table.read(path + 'THROUGHPUT/mse_itc_throughput_spec_mr.dat', format='ascii')
        elif spec.spectro == 'HR':
            spectro = Table.read(path + 'THROUGHPUT/mse_itc_throughput_spec_hr.dat', format='ascii')
            fiber = Table.read(path + 'THROUGHPUT/mse_itc_throughput_fits_hr.dat', format='ascii')

        # Interpolate
        thr_fiber = np.interp(spec.wgrid, fiber['lamA'], fiber['thr'])
        thr_poss = np.interp(spec.wgrid, poss['lamA'], poss['thr'])
        # Interpolate each arm for the spectrograph
        thr_spectro = np.zeros_like(thr_fiber)
        for i in range(1+int(np.max(spec.armgrid))):
            thr_spectro[spec.armgrid == i] = np.interp(spec.wgrid[spec.armgrid == i], spectro['lamA'],
                                                       spectro['thr_arm'+str(i+1)])

        # Apply throughput to sky and target
        spec.tgtflux *= thr_poss * thr_fiber * thr_spectro
        spec.skyflux *= thr_poss * thr_fiber * thr_spectro

        # Save throughput for verification
        spec.thr_poss = thr_poss
        spec.thr_fiber = thr_fiber
        spec.thr_spectro = thr_spectro

        return spec


    def compute_snr(self, doplot='online'):
        """Compute the SNR given a target spectrum and an observing configuration.

        Parameters
        ----------

        Returns
        -------
        """

        # Make a copy of the input spectrum so it's not modified
        s0 = copy.deepcopy(self)

        # Apply extinction
        s1 = s0.apply_atmos_ext()
        # Apply frontend throughput
        s2 = s1.apply_throughput_front()
        # Apply injection efficiency
        s3 = s2.apply_injeff()
        # Apply backend throughput
        s4 = s3.apply_throughput_back()

        # Surface of the telescope (60 segments of 1.44m corner to corner)
        surf = 60. * 3. / 2. * np.sqrt(3.) * (1.44 * 100. / 2.) ** 2.

        # Counts per second and per resolution element
        s4.tgtdetec = s4.tgtflux * surf * s4.reselgrid / s4.egrid
        s4.skydetec = s4.skyflux * surf * s4.reselgrid / s4.egrid

        # Resolution element
        if self.spectro == 'LR' or self.spectro == 'MR':
            npixspat = 4.1
            npixspec = 3.5
        elif self.spectro == 'HR':
            npixspat = 4.5
            npixspec = 4.5

        # Detector/spectrographs characteristics (per pixel, per second)
        if self.spectro == 'LR' or self.spectro == 'MR':
            dark = np.ones_like(s4.wgrid) * 0.02 / 3600.
            dark[s4.armgrid > 2] = 72. / 3600.
            readout = np.ones_like(s4.wgrid) * 5.
            readout[s4.armgrid > 2] = 8.
            nreadout = np.ones_like(s4.wgrid)
            if self.spectro == 'LR':
                nreadout[s4.armgrid > 0] = 2.  # (ASSUMING 2 READOUT IN GREEN ARM)
                nreadout[s4.armgrid > 1] = 4.  # (ASSUMING 4 READOUT IN RED ARM)
                nreadout[s4.armgrid > 2] = 12.  # (ASSUMING 12 READOUT UP THE RAMP IN NIR ARM)
            well = np.ones_like(s4.wgrid) * 70000.
            well[s4.armgrid > 2] = 45000.
            thermal = np.zeros_like(s4.wgrid)
            thermal[s4.armgrid > 2] = 9. / 3600.
        elif self.spectro == 'HR':
            dark = np.ones_like(s4.wgrid) * 0.02 / 3600.
            readout = np.ones_like(s4.wgrid) * 5.
            nreadout = np.ones_like(s4.wgrid)
            well = np.ones_like(s4.wgrid) * 70000.
            thermal = np.zeros_like(s4.wgrid)

        # Contamination from other source
        # -- need to create contaminating target with same conditions
        if self.badtgtmag > 0:
            badtgt0 = MseSpectrum(sessionID=-1, tgtmag=self.badtgtmag, band=self.band, template='flat', redshift=0,
                                  airmass=self.airmass, skymag=self.skymag, seeing=self.seeing, coating=self.coating,
                                  fibdiam=self.fibdiam, spectro=self.spectro, src_type='point')
            badtgt1 = badtgt0.apply_atmos_ext()
            badtgt2 = badtgt1.apply_throughput_front()
            badtgt3 = badtgt2.apply_injeff()
            badtgt4 = badtgt3.apply_throughput_back()
            badspec = badtgt4.tgtflux * surf * s4.reselgrid / s4.egrid

        else:
            badspec = np.zeros_like(s4.tgtdetec)  # best case is, sky is the only contaminant so bad-target is null
        # -- compute x-talk (bad target + sky)
        s4.xtalk = 0.001 * (badspec + s4.skydetec)
        # -- compute ghost (bad target + sky)
        s4.ghost = 0.001 * (badspec + s4.skydetec)
        # -- compute instrument diffuse light (real target + sky, wrong wavelength)
        s4.instdiffuse = np.zeros_like(s4.tgtdetec)
        for i in range(1 + int(np.max(s4.armgrid))):  # max within the same arm
            s4.instdiffuse[s4.armgrid == i] = 0.01 * np.max((s4.tgtdetec[s4.armgrid == i] + signal.medfilt(s4.skydetec[s4.armgrid == i], 101)))
        # -- compute telescope diffuse light (sky)
        s4.teldiffuse = np.zeros_like(s4.tgtdetec)
        for i in range(1 + int(np.max(s4.armgrid))):  # max within the same arm
            s4.teldiffuse[s4.armgrid == i] = 0.01 * signal.medfilt(s4.skydetec[s4.armgrid == i], 101)

        # Dark current per resolution element
        s4.dark = dark * npixspat * npixspec
        # Readout noise per resolution element
        s4.readout = readout * np.sqrt(npixspat / s4.spatbin * npixspec / s4.specbin * nreadout)
        # Thermal per resolution element
        s4.thermal = thermal * npixspat * npixspec

        # Compute SNR or exptime
        if s4.meth == 'getSNR':
            # account for exposure time
            s4.tgtdetec *= s4.etime
            s4.skydetec *= s4.etime
            s4.dark *= s4.etime
            s4.thermal *= s4.etime
            s4.xtalk *= s4.etime
            s4.ghost *= s4.etime
            s4.teldiffuse *= s4.etime
            s4.instdiffuse *= s4.etime
            # Photon noise
            s4.skynoise = np.sqrt(s4.skydetec)
            s4.tgtnoise = np.sqrt(s4.tgtdetec)
            s4.darknoise = np.sqrt(s4.dark)
            s4.thermalnoise = np.sqrt(s4.thermal)
            s4.xtalknoise = np.sqrt(s4.xtalk)
            s4.ghostnoise = np.sqrt(s4.ghost)
            s4.teldiffusenoise = np.sqrt(s4.teldiffuse)
            s4.instdiffusenoise = np.sqrt(s4.instdiffuse)
            # compute SNR
            s4.snr = s4.tgtdetec / np.sqrt(
                s4.tgtnoise ** 2 + s4.skynoise ** 2 + s4.darknoise ** 2 + s4.thermalnoise ** 2 + s4.readout ** 2 +
                s4.xtalknoise ** 2 + s4.ghostnoise ** 2 + s4.teldiffusenoise ** 2 + s4.instdiffusenoise ** 2)
        elif s4.meth == 'getEtime':
            aa = s4.tgtdetec ** 2
            bb = - s4.snr ** 2 * (s4.tgtdetec + s4.skydetec + s4.dark + s4.thermal + s4.xtalk + s4.ghost +
                                  s4.teldiffuse + s4.instdiffuse)
            cc = - s4.snr ** 2 * s4.readout
            s4.etime = (- bb + np.sqrt(bb ** 2 - 4. * aa * cc)) / (2. * aa)
            # account for exposure time
            s4.tgtdetec *= s4.etime
            s4.skydetec *= s4.etime
            s4.dark *= s4.etime
            s4.thermal *= s4.etime
            s4.xtalk *= s4.etime
            s4.ghost *= s4.etime
            s4.teldiffuse *= s4.etime
            s4.instdiffuse *= s4.etime
            # Photon noise
            s4.skynoise = np.sqrt(s4.skydetec)
            s4.tgtnoise = np.sqrt(s4.tgtdetec)
            s4.darknoise = np.sqrt(s4.dark)
            s4.thermalnoise = np.sqrt(s4.thermal)
            s4.xtalknoise = np.sqrt(s4.xtalk)
            s4.ghostnoise = np.sqrt(s4.ghost)
            s4.teldiffusenoise = np.sqrt(s4.teldiffuse)
            s4.instdiffusenoise = np.sqrt(s4.instdiffuse)

        # Plots
        output_file("plots.html")

        # Prepare figure
        if s4.meth == 'getSNR':
            fig1 = figure(title="SNR", x_axis_label="Wavelength (A)")
        else:
            fig1 = figure(title="Time", y_axis_type="log", x_axis_label="Wavelength (A)", y_axis_label="Seconds",
                          y_range=(np.nanmin(s4.etime)/1.1, np.nanmax(s4.etime)*1.1))
            etime_hours = s4.etime / 3600.
            fig1.extra_y_ranges = {"hours": Range1d(start=np.nanmin(etime_hours)/1.1, end=np.nanmax(etime_hours)*1.1)}
            fig1.add_layout(LogAxis(y_range_name="hours", axis_label="Hours"), 'right')

        if s4.meth == 'getSNR':
            # SNR and plot
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig1.line(s4.wgrid[arm], s4.snr[arm], line_color='black', line_alpha=.25)
                fig1.line(s4.wgrid[arm], signal.medfilt(s4.snr[arm], 101), line_color='black')
            if s4.spectro == 'LR':
                fig1.line([3600, 4000, 4000, 18000], [1, 1, 2, 2], line_color='cyan', line_dash="dashed")
            elif s4.spectro == 'MR':
                fig1.line([3600, 4000, 4000, 9500], [1, 1, 2, 2], line_color='cyan', line_dash="dashed")
            elif s4.spectro == 'HR':
                fig1.line([3600, 4000, 4000, 9000], [5, 5, 10, 10], line_color='cyan', line_dash="dashed")
        else:
            # Time and plot
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig1.line(s4.wgrid[arm], s4.etime[arm], line_color='black', line_alpha=.25)
                fig1.line(s4.wgrid[arm], signal.medfilt(s4.etime[arm], 101), line_color='black')

        if doplot == 'online':
            script, div = components(fig1)

        elif doplot == 'offline':
            fig2 = figure(title="Spectra", y_axis_type="log", y_axis_label="Flux (erg/s/cm2/A)", x_axis_label="Wavelength (A)")
            fig3 = figure(title="Counts", y_axis_type="log", y_axis_label="Counts (photons/s/cm2/A)", x_axis_label="Wavelength (A)")
            fig4 = figure(title="Throughput", x_axis_label="Wavelength (A)")
            fig6 = figure(title="Noise", y_axis_type="log", y_axis_label="Counts (photons/s/cm2/A)", x_axis_label="Wavelength (A)")

            # Plot intrinsic spectra
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig2.line(s0.wgrid[arm], s0.tgtflux[arm], line_color='#FFBB00', line_alpha=.25)
                fig2.line(s0.wgrid[arm], s0.skyflux[arm], line_color='#0088BB', line_alpha=.25)
                # Overline spectrum after extinction
                fig2.line(s1.wgrid[arm], s1.tgtflux[arm], line_color='#DD8800', line_alpha=.25)
                # Overline spectrum after throughput+injection
                fig2.line(s4.wgrid[arm], s4.tgtflux[arm], line_color='#FF0000', line_alpha=.25)
                fig2.line(s4.wgrid[arm], s4.skyflux[arm], line_color='#0000FF', line_alpha=.25)
                # Overline median filtered spectra
                fig2.line(s0.wgrid[arm], signal.medfilt(s0.tgtflux[arm], 101), line_color='#FFBB00', legend='Target')
                fig2.line(s0.wgrid[arm], signal.medfilt(s0.skyflux[arm], 101), line_color='#0088BB', legend='Sky')
                fig2.line(s1.wgrid[arm], signal.medfilt(s1.tgtflux[arm], 101), line_color='#DD8800', legend='Target + atmosphere')
                fig2.line(s4.wgrid[arm], signal.medfilt(s4.tgtflux[arm], 101), line_color='#FF0000', legend='Target out')
                fig2.line(s4.wgrid[arm], signal.medfilt(s4.skyflux[arm], 101), line_color='#0000FF', legend='Sky out')

            # Plot counts on detector
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig3.line(s4.wgrid[arm], s4.dark[arm], line_color='#00FF00', legend='Dark')
                fig3.line(s4.wgrid[arm], s4.tgtdetec[arm], line_color='#FF0000', line_alpha=.25)
                fig3.line(s4.wgrid[arm], s4.skydetec[arm], line_color='#0000FF', line_alpha=.25)
                # Overline median filtered spectra
                fig3.line(s4.wgrid[arm], signal.medfilt(s4.tgtdetec[arm], 101), line_color='#FF0000', legend='Target counts')
                fig3.line(s4.wgrid[arm], signal.medfilt(s4.skydetec[arm], 101), line_color='#0000FF', legend='Sky counts')

            # Throughput plot
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig4.line(s4.wgrid[arm], s4.thr_struc[arm], line_color='#FF0000', legend='Structure')
                fig4.line(s4.wgrid[arm], s4.thr_m1[arm], line_color='#0000FF', legend='M1')
                fig4.line(s4.wgrid[arm], s4.thr_pfue[arm], line_color='#AA4400', legend='PFUE')
                fig4.line(s4.wgrid[arm], s4.inj[arm], line_color='#00AA66', legend='Inj.Eff.')
                fig4.line(s4.wgrid[arm], s4.thr_poss[arm], line_color='#00FF88', legend='PosS')
                fig4.line(s4.wgrid[arm], s4.thr_fiber[arm], line_color='#8800FF', legend='FiTS')
                fig4.line(s4.wgrid[arm], s4.thr_spectro[arm], line_color='#CCCC00', legend='Spectro')
                # overall throughput
                fig4.line(s4.wgrid[arm], (s4.thr_struc * s4.thr_m1 * s4.thr_pfue * s4.thr_poss * s4.inj * s4.thr_fiber
                                          * s4.thr_spectro)[arm], line_color='black')

            # Dark, readout, Poisson noise ...
            for i in range(1 + int(np.max(s4.armgrid))):
                arm = s4.armgrid == i
                fig6.line(s4.wgrid[arm], s4.darknoise[arm], line_color='#00FF00', legend='Dark noise')
                fig6.line(s4.wgrid[arm], s4.readout[arm], line_color='#FF8800', legend='Read noise')
                fig6.line(s4.wgrid[arm], s4.tgtnoise[arm], line_color='#FF0000', alpha=.25)
                fig6.line(s4.wgrid[arm], s4.skynoise[arm], line_color='#0000FF', alpha=.25)
                fig6.line(s4.wgrid[arm], signal.medfilt(s4.tgtnoise[arm], 101), line_color='#FF0000', legend='Target noise')
                fig6.line(s4.wgrid[arm], signal.medfilt(s4.skynoise[arm], 101), line_color='#0000FF', legend='Sky noise')
                fig6.line(s4.wgrid[arm], s4.thermalnoise[arm], line_color='#FF4488', legend='Thermal noise')
                fig6.line(s4.wgrid[arm], s4.xtalknoise[arm], line_color='#004488', legend='X-talk noise')
                fig6.line(s4.wgrid[arm], s4.instdiffusenoise[arm], line_color='#440088', legend='Diffuse (inst.) noise')
                fig6.line(s4.wgrid[arm], s4.teldiffusenoise[arm], line_color='#00FF88', legend='Diffuse (tel.) noise')

            # make a grid
            grid = gridplot([[fig1, fig2, fig6], [fig4, None, fig3]])
            show(grid)
            script, div = components(grid)

        else:
            script = ''
            div = ''

        return s4, script, div

def mag2flux(lamA, mag):
    """Returns the flux and counts at a given wavelength for a given AB magnitude.
    """

    # flux in Jansky: 1Jy=10-23 ergs/s/cm2/Hz
    fnu = 10 ** (-(mag - 8.9) / 2.5)
    # flux in ergs/s/cm2/A: -23-10+20=-13
    flam = fnu * const.c.value * 1e-13 / lamA ** 2.
    # energy of photons, in ergs (1J=10^7ergs, 1A=10^-10m)
    phot_e = 1e17 * const.h.value * const.c.value / lamA
    # number of photons /s/cm2/A
    count = flam / phot_e

    return flam, count


def flux2mag(lamA, flux, filtrans):
    """Returns the AB magnitude for a given filter's transmission and a given spectrum.

    Parameters
    ----------
    flux -- the spectrum flux in ergs/s/cm2/A
    wave -- the spectrum wavelength in A
    filtrans -- the filter transmission (1 for 100%)

    Returns
    -------
    mag -- the AB magnitude
    """
    # flux from ergs/s/cm2/A to Jy: -23-10+20=-13
    flux_tmp = flux / (const.c.value * 1e-13 / lamA ** 2.)
    # energy of photons, in ergs (1J=10^7ergs, 1A=10^-10m)
    phot_e = 1e17 * const.h.value * const.c.value / lamA
    mag = -2.5 * np.log10(np.trapz(flux_tmp / phot_e * filtrans, x=lamA) / np.trapz(3631. / phot_e * filtrans, x=lamA))

    return mag

