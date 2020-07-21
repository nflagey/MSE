#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, vstack
import copy
from bokeh.plotting import figure, save
from astropy import wcs

# Global variables
patrol = 1.24 * 7.77 / 106.7e-3  # in arcsec (average plate scale)


class MseFibAlloc:

    def __init__(self, file='targets.dat', spectro='HR', meth='one', iternum=5, allocfrac=95,
                 fovctr='auto', fovctrra='180', fovctrdec='180', doplot=False, dither=False):

        # Input file
        self.file = file

        # Making plots? (very time consuming)
        self.doplot = doplot

        # Dithering?
        self.dither = dither

        # Input parameters
        self.spectro = spectro
        self.meth = meth
        self.iternum = iternum
        self.allocfrac = allocfrac
        self.fovctr = fovctr
        self.fovctrra = fovctrra
        self.fovctrdec = fovctrdec

        # Targets info
        self.tgt_in = Table()  # input
        self.tgt_fov = Table()  # can be reached
        self.tgt_wip = Table()  # working copy

        # Fibers info
        self.fib_in = Table()
        self.fib_fov = Table()
        self.fib_wip = Table()
        # Distances
        self.dist = []      # this is for all fibers and all targets
        self.dist_fov = []  # this is for all fibers and all targets within reach of at least 1 fiber
        self.dist_wip = []  # this is for ...
        self.all_dist = {}  # dictionary with all distances for each iteration
        self.pairs = []
        self.trackfrac = []
        self.kmax = 10

        # Create targets
        self.create_tgt()
        # Create fibers
        self.create_fib()
        # Compute distances
        self.compute_dist()
        # Run the optimization
        self.run_optim()

    def create_tgt(self):
        # Open CSV file with Targets (input file should have the following columns:
        # RAJ2000, DECJ2000, umag, gmag, rmag, imag, zmag, Jmag, Hmag, Nobsreq, Nrepeat, priority, surveypriority)
        tgt = Table.read('TARGETS/' + self.file, format='csv')  # TODO: accept multiple input files

        # Project targets coordinates onto field of view
        im = np.empty([360, 360])
        hdu = fits.PrimaryHDU(im)
        hdul = fits.HDUList([hdu])
        # Create default header
        hdr = hdul[0].header
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRPIX1'] = 0.
        hdr['CRPIX2'] = 0.
        hdr['CDELT1'] = -1. / 3600
        hdr['CDELT2'] = 1. / 3600
        # Add coordinates (depends on input file)
        if self.fovctr == 'manual':
            hdr['CRVAL1'] = self.fovctrra
            hdr['CRVAL2'] = self.fovctrdec
        else:
            hdr['CRVAL1'] = np.mean(tgt['RAJ2000'])
            hdr['CRVAL2'] = np.mean(tgt['DECJ2000'])
        # Get pixel coordinates for all targets
        wcs0 = wcs.WCS(hdul[0])
        xy = [wcs0.wcs_world2pix([[tgt['RAJ2000'][i], tgt['DECJ2000'][i]]], 1)
              for i in range(len(tgt['RAJ2000']))]
        # Make those columns
        xpos = Column([xy[i][0][0] for i in range(len(xy))], name='xpos', unit='arcsec')
        ypos = Column([xy[i][0][1] for i in range(len(xy))], name='ypos', unit='arcsec')
        # Add an "Nobsdone" column if it does not exist already
        if 'Nobsdone' not in tgt.colnames:
            nobs = Column([0 for i in range(len(xy))], name='Nobsdone')
            # Add columns to table
            tgt.add_columns([xpos, ypos, nobs])
        else:
            tgt.add_columns([xpos, ypos])
        # Store
        self.tgt_in = tgt

    def create_fib(self):
        # read AAO file (coordinates in mm)
        fib = Table.read('POSS/' + 'sphinx.tsv', format='ascii')
        # convert into arcsec
        fib['col1'] /= 106.7e-3
        fib['col2'] /= 106.7e-3

        # Select PosS type and create a Table (X, Y, Dither)
        fib_spec = Table([fib[fib['col3'] == self.spectro]['col1'],
                          fib[fib['col3'] == self.spectro]['col2'],
                          fib[fib['col3'] == self.spectro]['col2'] * 0],  # dithering position 0
                         names=('xpos', 'ypos', 'dither'))
        # Dithering?
        if self.dither:
            print("<h3>You are using Dithering<h3>")
            fib_spec = vstack([fib_spec, Table([fib[fib['col3'] == self.spectro]['col1'] - 7.77/106.7e-3,
                                                fib[fib['col3'] == self.spectro]['col2'] + 7.77*np.sqrt(3)/6/106.7e-3,
                                                fib[fib['col3'] == self.spectro]['col2'] * 0 + 1],  # dith pos 1
                                               names=('xpos', 'ypos', 'dither'))])
            fib_spec = vstack([fib_spec, Table([fib[fib['col3'] == self.spectro]['col1'],
                                                fib[fib['col3'] == self.spectro]['col2'] + 7.77*np.sqrt(3)/3/106.7e-3,
                                                fib[fib['col3'] == self.spectro]['col2'] * 0 + 2],  # dith pos 2
                                               names=('xpos', 'ypos', 'dither'))])

        self.fib_in = fib_spec

    def compute_dist(self):

        # Initiate distances array
        ntgt_in = len(self.tgt_in)
        nfib_in = len(self.fib_in)
        # Put all distances at 666 to begin with
        distances = np.ones((ntgt_in, nfib_in)) * 666

        # Loop over fibers
        for f in range(nfib_in):
            # only select targets that are within a small box (+/-10mm) around the fiber
            # since the patrol region is 9.64mm radius
            near = ((np.abs(self.tgt_in['xpos'] - self.fib_in['xpos'][f]) <= 10 / 106.7e-3)
                    & (np.abs(self.tgt_in['ypos'] - self.fib_in['ypos'][f]) <= 10 / 106.7e-3))
            targets_near = self.tgt_in[near]
            lnear = len(targets_near)
            # compute distances
            if lnear > 0:
                distances[near, f] = [np.sqrt((targets_near['xpos'][i] - self.fib_in['xpos'][f]) ** 2
                                              + (targets_near['ypos'][i] - self.fib_in['ypos'][f]) ** 2)
                                      for i in range(lnear)]

        # Put to 666 all distances larger than patrol region
        distances[distances > patrol] = 666
        # Store distances
        self.dist = distances

        # Remove targets that cannot be reached by any fiber and fiber that cannot reach any fiber
        # -- targets
        sel_tgt_fov = np.min(self.dist, axis=1) != 666
        self.tgt_fov = self.tgt_in[sel_tgt_fov]
        # add "FoV" column (True/False) to check results
        self.tgt_in.add_column(Column(sel_tgt_fov, dtype=int), name="fov")
        # -- fibers
        sel_fib_fov = np.min(self.dist, axis=0) != 666
        self.fib_fov = self.fib_in[sel_fib_fov]
        # -- distances
        self.dist_fov = self.dist[sel_tgt_fov, :][:, sel_fib_fov]

        # Number of fibers/targets to work with
        ntgt_fov = len(self.tgt_fov)

        print(" <h2> Initial conditions </h2> ")
        print(str(ntgt_fov) + " targets can be reached by a positioner out of " + str(
            ntgt_in) + " targets provided by user</br></br>")
        print(" <h2> Optimizing </h2> ")

    def run_optim(self):
        # Run optimization
        self.niter = 0
        trackfrac = []  # save fraction of allocated targets within FoV
        while True:
            print('<h3> Iteration #'+str(self.niter+1)+'</h3')
            # Run one allocation
            self.optim()
            # Compute fraction of targets observed so far
            trackfrac = np.append(trackfrac, 1. * len(self.tgt_fov[self.tgt_fov['Nrepeat'] == self.tgt_fov['Nobsdone']])
                                  / len(self.tgt_fov))
            # Print how many targets are left
            print(str(len(self.tgt_fov[self.tgt_fov['Nrepeat'] != self.tgt_fov['Nobsdone']])) + " targets left out of "
                  + str(len(self.tgt_fov)) + " targets to observe </br>")
            # Increased number of iterations
            self.niter += 1
            # If all possible targets have been observed, then stop
            if 1. * len(self.tgt_fov[self.tgt_fov['Nrepeat'] == self.tgt_fov['Nobsdone']]) == len(self.tgt_fov):
                break
            # What is the computing method?
            if self.meth == "one":
                if self.niter == 1:
                    break
            elif self.meth == "fixiter":
                if self.niter == self.iternum:
                    break
            elif self.meth == "fixgoal":
                if trackfrac[-1] >= self.allocfrac/100.:
                    break

        # Update Target table and print into file
        sel_tgt_fov = np.min(self.dist, axis=1) != 666
        self.tgt_in[sel_tgt_fov] = self.tgt_fov
        self.tgt_in.write('TARGETS/' + 'results.csv', format='csv', overwrite=True)

        self.trackfrac = trackfrac

    def optim(self):

        # 0. Select correct FIB/DIST because of DITHERING
        if self.dither:
            fib_fov = self.fib_fov[self.fib_fov['dither'] == (self.niter % 3)]
            dist_fov = self.dist_fov[:, self.fib_fov['dither'] == (self.niter % 3)]
        else:
            fib_fov = self.fib_fov[self.fib_fov['dither'] == 0]
            dist_fov = self.dist_fov[:, self.fib_fov['dither'] == 0]

        # 1. Build a working copy of the important variables to keep the original un-modified

        # Number of fibers/targets (FoV) once for all
        ntgt_fov = len(self.tgt_fov)
        nfib_fov = len(fib_fov)

        # Remove targets that cannot be reached by any fiber and fiber that cannot reach any fiber
        sel_tgt_wip = (np.min(dist_fov, axis=1) != 666) & (self.tgt_fov['Nrepeat'] != self.tgt_fov['Nobsdone'])
        # Apply selection for targets first otherwise filtering fibers won't work
        self.tgt_wip = self.tgt_fov[sel_tgt_wip]

        # Check that there are actually targets that can be reached (need to do that test here for the DITHERING case)
        if len(self.tgt_wip) == 0:
            print("No target can be reached by any positioner in this dithering position but there are targets left."
                  "Going to next dithering position </br>")
            # Keep all allocated distances in memory (need to keep that for the plots I make outside)
            self.all_dist[self.niter] = self.fib_wip['dist'][0] * 0 - 666
        else:
            self.dist_wip = dist_fov[sel_tgt_wip, :]
            # Remove fibers that cannot reach any fiber
            sel_fib_wip = np.min(self.dist_wip, axis=0) != 666
            # Apply selection
            self.fib_wip = fib_fov[sel_fib_wip]
            self.dist_wip = self.dist_wip[:, sel_fib_wip]
            # Number of fibers/targets to work with
            ntgt_wip = len(self.tgt_wip)
            nfib_wip = len(self.fib_wip)

            print(str(ntgt_wip) + " targets can be reached by a positioner out of " + str(ntgt_fov)
                  + " targets to observe </br>")
            print(str(nfib_wip) + " fibers can reach a target out of " + str(nfib_fov) + " fibers </br></br>")

            # 2. Create LONG/SHORT variables so to switch easily between fibers and targets

            # Transpose if ntgt < nfib
            if ntgt_wip < nfib_wip:
                # Create the relevant variables
                short, nshort, wshort = self.tgt_wip, ntgt_wip, 'targets'
                long, nlong, wlong = self.fib_wip, nfib_wip, 'fibers'
                # Working copy of distance
                dist_wip = np.transpose(self.dist_wip)
                # Apply priority
                for t in range(nshort):
                    dist_wip[:, t] = self.dist_wip[t, :] / self.tgt_wip['priority'][t]
                dist_wip[np.transpose(self.dist_wip) == 666] = 666  # make sure targets outside patrol region remain outside
            else:
                short, nshort, wshort = self.fib_wip, nfib_wip, 'fibers'
                long, nlong, wlong = self.tgt_wip, ntgt_wip, 'targets'
                # Working copy of distance
                dist_wip = copy.deepcopy(self.dist_wip)
                # Apply priority
                for t in range(nlong):
                    dist_wip[t, :] = self.dist_wip[t, :] / self.tgt_wip['priority'][t]
                dist_wip[self.dist_wip == 666] = 666  # make sure those targets outside patrol region remain outside

            # Matching array
            state = np.zeros_like(dist_wip)
            # Find nearest target to each fiber
            best = np.argmin(dist_wip, axis=0)
            # Create matching array
            state[best, range(nshort)] = 1
            # Get pairs
            test = [(0 < i) & (i < patrol) for i in state * dist_wip]
            pairs = list(np.where(test))

            # Any two fibers paired to the same target? (or two targets to the same fiber)
            short_per_long = np.sum(state, axis=1)
            more_than_one_short = list(filter(lambda x: x > 1, short_per_long))
            # --> yes, some doubles ...
            if len(more_than_one_short) > 0:

                # Set up initial energy and store, at least for tests
                energy = 666. * nshort
                all_energy = np.array([energy])

                # Loop a few times on the simulation
                for k in range(self.kmax):

                    # Working copy of STATE and DIST
                    state_wip = np.zeros_like(dist_wip)
                    dist_wip2 = copy.deepcopy(dist_wip)

                    # define random order
                    tmp = np.argsort(np.random.random_sample(nshort))
                    # loop over fibers in that order
                    for f in tmp:
                        nearest = np.argmin(dist_wip2[:, f])
                        state_wip[nearest, f] = 1
                        # put that target's distance at 666 for all other fibers so it won't be allocated again
                        sav = dist_wip2[nearest, f]
                        dist_wip2[nearest, :] = 666
                        dist_wip2[nearest, f] = sav

                    # Compute energy for that state
                    energy_wip = np.sum(dist_wip2 * state_wip)

                    # Pairs
                    test = [(0 < i) & (i < patrol) for i in state_wip * dist_wip2]
                    pairs_wip = np.asarray(np.where(test))

                    # Compute weight for wide/narrow range of magnitudes
                    if ntgt_wip < nfib_wip:
                        minumag, maxumag = np.nanmin(short['umag'][pairs_wip[0]]), np.nanmax(short['umag'][pairs_wip[0]])
                        mingmag, maxgmag = np.nanmin(short['gmag'][pairs_wip[0]]), np.nanmax(short['gmag'][pairs_wip[0]])
                        minrmag, maxrmag = np.nanmin(short['rmag'][pairs_wip[0]]), np.nanmax(short['rmag'][pairs_wip[0]])
                        minimag, maximag = np.nanmin(short['imag'][pairs_wip[0]]), np.nanmax(short['imag'][pairs_wip[0]])
                        minzmag, maxzmag = np.nanmin(short['zmag'][pairs_wip[0]]), np.nanmax(short['zmag'][pairs_wip[0]])
                        minjmag, maxjmag = np.nanmin(short['Jmag'][pairs_wip[0]]), np.nanmax(short['Jmag'][pairs_wip[0]])
                        minhmag, maxhmag = np.nanmin(short['Hmag'][pairs_wip[0]]), np.nanmax(short['Hmag'][pairs_wip[0]])
                    else:
                        minumag, maxumag = np.nanmin(long['umag'][pairs_wip[1]]), np.nanmax(long['umag'][pairs_wip[1]])
                        mingmag, maxgmag = np.nanmin(long['gmag'][pairs_wip[1]]), np.nanmax(long['gmag'][pairs_wip[1]])
                        minrmag, maxrmag = np.nanmin(long['rmag'][pairs_wip[1]]), np.nanmax(long['rmag'][pairs_wip[1]])
                        minimag, maximag = np.nanmin(long['imag'][pairs_wip[1]]), np.nanmax(long['imag'][pairs_wip[1]])
                        minzmag, maxzmag = np.nanmin(long['zmag'][pairs_wip[1]]), np.nanmax(long['zmag'][pairs_wip[1]])
                        minjmag, maxjmag = np.nanmin(long['Jmag'][pairs_wip[1]]), np.nanmax(long['Jmag'][pairs_wip[1]])
                        minhmag, maxhmag = np.nanmin(long['Hmag'][pairs_wip[1]]), np.nanmax(long['Hmag'][pairs_wip[1]])
                    # Average of magnitude range over 7 bands (+ 1 to avoid getting 0)
                    magscale = ((maxumag - minumag) + (maxgmag - mingmag) + (maxrmag - minrmag) + (maximag - minimag)
                                + (maxzmag - minzmag) + (maxjmag - minjmag) + (maxhmag - minhmag)) / 7. + 1.

                    print(magscale)

                    # Store current value
                    all_energy = np.append(all_energy, energy_wip / magscale)

                    # Is it the lowest energy?
                    if energy_wip < energy:
                        print("keeper", magscale)
                        energy = copy.deepcopy(energy_wip)
                        pairs = copy.deepcopy(pairs_wip)

                # End loop on the simulation

            # Final pairs
            npairs = len(pairs[0])
            print(str(npairs) + " pairs out of " + str(nshort) + " " + wshort +
                  " and " + str(nlong) + " " + wlong + ".</br>")

            # leaving short/long space and going back to fiber/target space
            if nfib_wip > ntgt_wip:
                pairs[0], pairs[1] = pairs[1], pairs[0]
            # find what targets/fibers are paired and store in FIB table
            self.fib_wip['xtgt'] = [None] * nfib_wip
            self.fib_wip['ytgt'] = [None] * nfib_wip
            self.fib_wip['dist'] = [666.] * nfib_wip
            self.pairs = pairs
            for i in range(npairs):
                self.fib_wip['xtgt'][pairs[1][i]] = self.tgt_wip['xpos'][pairs[0][i]]
                self.fib_wip['ytgt'][pairs[1][i]] = self.tgt_wip['ypos'][pairs[0][i]]
                self.fib_wip['dist'][pairs[1][i]] = self.dist_wip[pairs[0][i], pairs[1][i]]  # store real distances

            # Add 1 to NOBS
            self.tgt_wip['Nobsdone'][pairs[0]] += 1
            self.tgt_fov[sel_tgt_wip] = self.tgt_wip

            # Print fraction of targets observed so far
            frac = 100. * len(self.tgt_fov[self.tgt_fov['Nrepeat'] == self.tgt_fov['Nobsdone']]) / len(self.tgt_fov)
            print('Allocation fraction so far: ' + '{:.1f}'.format(frac) + '%</br></br>')

            # Keep all allocated distances in memory
            self.all_dist[self.niter] = self.fib_wip['dist'][pairs[1]]

            # Plot field of view (really slow on entire field of view) !!!
            if self.doplot:
                fig = figure(title="Final allocation map (" + self.spectro + ")",
                             x_axis_label="RA offset (arcsec)", y_axis_label="DEC offset (arcsec)",
                             width=800, height=800, tools='pan, wheel_zoom, zoom_in, zoom_out, reset, save',
                             match_aspect=True)
                self.optim_plot(fig)
                save(fig)

    def optim_plot(self, fig):

        # Plot field of view
        npairs = len(self.pairs[0])
        # All PosS
        if self.spectro == 'HR':
            color = 'mediumpurple'
        else:
            color = 'crimson'
        fig.circle(self.fib_in['xpos'], self.fib_in['ypos'],
                   color=color, alpha=0.5, radius=patrol)
        # All targets provided
        fig.x(self.tgt_in['xpos'], self.tgt_in['ypos'],
              color='black', alpha=0.5, legend='All targets')
        # Targets in FOV
        fig.x(self.tgt_fov['xpos'], self.tgt_fov['ypos'],
              color='black', line_width=2., legend='Targets not observed')
        # Targets observed
        fig.x(self.tgt_fov[self.tgt_fov['Nrepeat'] == self.tgt_fov['Nobsdone']]['xpos'],
              self.tgt_fov[self.tgt_fov['Nrepeat'] == self.tgt_fov['Nobsdone']]['ypos'],
              color='greenyellow', line_width=2, legend='Targets observed')

        # Allocated targets
        for i in range(npairs):
            fig.x(self.tgt_wip['xpos'][self.pairs[0][i]], self.tgt_wip['ypos'][self.pairs[0][i]],
                  color='red', size=6, line_width=2, legend='Targets allocated')
        # Pairs
        fig.multi_line([[self.fib_wip['xpos'][self.pairs[1][i]], self.tgt_wip['xpos'][self.pairs[0][i]]]
                        for i in range(npairs)],
                       [[self.fib_wip['ypos'][self.pairs[1][i]], self.tgt_wip['ypos'][self.pairs[0][i]]]
                        for i in range(npairs)],
                       color='red')

