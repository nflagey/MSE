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
import time
from scipy.spatial import cKDTree, KDTree
from sklearn.neighbors import BallTree
from scipy.sparse import lil_matrix, coo_matrix
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import random

# Global variables
plate_scale = 106.7e-3  # microns per arcsec
patrol = 1.24 * 7.77 / plate_scale  # in arcsec
bigR = 7.77 / plate_scale  # circumradius
smallR = np.sqrt(3) / 2 * bigR  # inradius

score_none = 666.  # distance for unallocated fiber's score

# GA parameters
init_pop_size = 20  # size of initial population
elite_pct = 25.  # percent of best parents to keep for next generation
mating_pct = 50.  # percent of best parents allowed to mate
mutation_pct = 5.  # percent of genes that mutate after crossover
max_generation = 10  # maximum number of generations


class MseFibAlloc(object):
    """
    A class that describes the allocation of fiber positioners to targets
    """

    def __init__(self, file='targets.dat', spectro='HR', meth='one', iternum=5, allocfrac=95,
                 fovctr='auto', fovctrra='180', fovctrdec='180', doplot=False, dither=False,
                 init_size=20, mat_size=10, off_size=10, pmutation=1, ngeneration=100):

        # SessionID
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

        # GA parameters
        self.init_size = init_size
        self.mat_size = mat_size
        self.off_size = off_size
        self.pmutation = pmutation
        self.ngeneration = ngeneration

        self._targets = None
        self._poss = None
        self._allocMatrix = None
        self._allocMatrix = None

    @property
    def poss(self):
        """
        Defines the fiber Positioners system
        :return: Table describing the Positioners system with four column (xpos, ypos, type, dither)
        """
        if getattr(self, '_poss', None) is None:
            # Read AAO file (coordinates in mm)
            poss = Table.read('POSS/' + 'sphinx.tsv', format='ascii')
            # only select LR or HR
            poss = poss[poss['col3'] == self.spectro]

            # only use a small central region
            #poss = poss[(poss['col1'] < 25) & (poss['col1'] > -25)
            #            & (poss['col2'] < 25) & (poss['col2'] > -25)]

            # convert into arcsec
            poss['col1'] /= plate_scale
            poss['col2'] /= plate_scale
            # rename columns
            poss['col1'].name = 'xpos'
            poss['col2'].name = 'ypos'
            poss['col3'].name = 'type'
            # add dither column
            poss['dither'] = 0
            print("<h3> You are using the {0} positioners.<h3><br>".format(self.spectro))
            # Dithering?
            if self.dither:
                print("<h3> You are using Dithering<h3>")
                poss = vstack([poss,
                               Table([poss['xpos'] + bigR / 4,
                                      poss['ypos'] + smallR / 2,
                                      poss['type'],
                                      poss['dither'] + 1],
                                     names=('xpos', 'ypos', 'type', 'dither')),
                               Table([poss['xpos'] + 3 * bigR / 4,
                                      poss['ypos'] + smallR / 2,
                                      poss['type'],
                                      poss['dither'] + 2],
                                     names=('xpos', 'ypos', 'type', 'dither')),
                               Table([poss['xpos'] + bigR / 2,
                                      poss['ypos'],
                                      poss['type'],
                                      poss['dither'] + 3],
                                     names=('xpos', 'ypos', 'type', 'dither'))])
            # pass it to self
            self._poss = poss
        return self._poss

    @property
    def targets(self):
        """
        Defines the Targets
        :return: Table
        """
        if getattr(self, '_targets', None) is None:

            # Read user's input file (need to define a format)
            tgt = Table.read('TARGETS/' + self.file, format='csv')
            #tgt['RAJ2000'] = tgt['RAJ2000'].astype(np.float32)
            #tgt['DECJ2000'] = tgt['DECJ2000'].astype(np.float32)
            #tgt['repeat'] = tgt['repeat'].astype(np.uint8)
            #tgt['priority'] = tgt['priority'].astype(np.uint8)

            # Project targets coordinates onto field of view
            im = np.empty([1, 1])
            hdu = fits.PrimaryHDU(im)
            hdul = fits.HDUList([hdu])
            # Create default header
            hdr = hdul[0].header
            hdr['CTYPE1'] = 'RA---TAN'
            hdr['CTYPE2'] = 'DEC--TAN'
            hdr['CRPIX1'] = 0.
            hdr['CRPIX2'] = 0.
            hdr['CDELT1'] = -1. / 3600  # 1 arcsec pixels
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
            nobs = Column(np.zeros(len(xy)), name='nobs')

            # Add columns to table
            tgt.add_columns([xpos, ypos, nobs])
            # Store
            self._targets = tgt

            # Print stats on screen
            print("<h3> Some statistics about the targets: <h3>")
            print("Number of targets entered by user: {0} ".format(len(tgt)))
            print("Number of observations (total(targets x repeat)) entered by user: {0} ".format(np.sum(tgt['Nrepeat'])))

        return self._targets

    @property
    def alloc_matrix_poss(self):
        """
        Find which Targets can be reached by each Positioner using KDTree
        :return:
        """
        if getattr(self, '_alloc_matrix_poss', None) is None:

            print("<h3> Identifying nearby targets ... <h3>")
            # extract XY columns only for each of the LR and HR
            poss_xy = np.array([self.poss['xpos'], self.poss['ypos']]).transpose()
            targets_xy = np.array([self.targets['xpos'], self.targets['ypos']]).transpose()
            # get kd-trees
            kdtTarget = cKDTree(targets_xy)  # fibers to targets
            # get nearby targets for each fiber
            indTarget = kdtTarget.query_ball_point(poss_xy, patrol, p=2)
            nTarget = np.array([len(indTarget[i]) for i in range(len(self.poss))])
            # --> this returns a list of lists (super fast!)

            # Compute distances and put everything in a nice Table
            distTarget = [np.array([np.sqrt((self.targets['xpos'][t] - self.poss['xpos'][f]) ** 2
                                            + (self.targets['ypos'][t] - self.poss['ypos'][f]) ** 2)
                                    for t in indTarget[f] if len(indTarget[f]) > 0]) for f in range(len(self.poss))]

            # populate matrix
            alloc_matrix_poss = Table([self.poss['xpos'], self.poss['ypos'], nTarget, indTarget, distTarget],
                                      names=['xpos', 'ypos', 'ntarget', 'reachable', 'distance'])
            # remove unreachable positioners
            alloc_matrix_poss = alloc_matrix_poss[nTarget > 0]
            # pass to self
            self._alloc_matrix_poss = alloc_matrix_poss

            # alloc = [np.zeros(len(indTarget[f]), dtype='bool') for f in range(len(poss))]
            # alloc_1 = [None for f in range(len(poss))]  # this gives the target index in the target list
            # alloc_2 = [None for f in range(len(poss))]  # this gives the target index in the reachable list
            # prior = [[self.targets['priority'][t] for t in indTarget[f]] for f in range(len(poss))]
            # repeat = [[self.targets['repeat'][t] for t in indTarget[f]] for f in range(len(poss))]
            # nobs = [np.zeros(len(indTarget[f]), dtype='int') for f in range(len(poss))]
            # dith = self.poss['dither']

        return self._alloc_matrix_poss

    @property
    def alloc_matrix_target(self):
        """
        Find which Positioner can be allocated to each Target using KDTree
        :return:
        """
        if getattr(self, '_alloc_matrix_target', None) is None:
            print("<h3> Identifying nearby targets ... <h3>")
            # extract XY columns only for each of the LR and HR
            poss = np.array([self.poss['xpos'], self.poss['ypos']]).transpose()
            targets = np.array([self.targets['xpos'], self.targets['ypos']]).transpose()
            # get kd-trees
            kdtPoss = cKDTree(poss)  # targets to fibers
            # get nearby targets for each fiber
            indPoss = kdtPoss.query_ball_point(targets, patrol, p=2)
            # --> this returns a list of lists (super fast!)

            # Compute distances and put everything in a nice Table
            distPoss = [np.array([np.sqrt((self.poss['xpos'][f] - self.targets['xpos'][t]) ** 2
                                          + (self.poss['ypos'][f] - self.targets['ypos'][t]) ** 2)
                                  for f in indPoss[t] if len(indPoss[t]) > 0]) for t in range(len(targets))]
            # dith = self.poss['dither']

            # populate matrix
            alloc_matrix_target = Table([indPoss, distPoss], names=['reachable', 'distance'])
            # pass to self
            self._alloc_matrix_target = alloc_matrix_target

        return self._alloc_matrix_target

    def ga_loop(self):
        """
        Generate the initial population by allocating randomly targets to fibers
        :return:
        """

        t0 = time.time()

        # === Initial population ===
        # only pick solution that are allowed, using alloc_matrix_poss
        init_pop = np.array([np.random.choice(self.alloc_matrix_poss['reachable'][i], self.init_size)
                             for i in range(len(self.alloc_matrix_poss))])
        # this is an np.array (nposs, pop_size)
        print(init_pop)

        t1 = time.time()
        print(f"Initial population (ms): {1e3*(t1-t0):.3f}")

        # === Evaluate initial population ===
        init_scores = self.ga_evaluate(init_pop)
        # this is a np.array (nposs, pop_size)
        init_score = np.sum(init_scores, axis=0)
        # this is a np.array (pop_size)

        # save best score
        self.nalloc = np.zeros(self.ngeneration + 1)
        self.best = np.zeros(self.ngeneration + 1)
        self.best[0] = np.min(init_score)
        # find best chromosome for plot
        self.solution = init_pop[:, np.argmin(init_score)]
        best_scores = init_scores[:, np.argmin(init_score)]
        self.nalloc[0] = len(best_scores[best_scores != 666])

        t2 = time.time()
        print(f"Evaluate initial population (ms): {1e3*(t2 - t1): .3f}")

        # === Loop to create new generations ===
        for i in range(self.ngeneration):
            ta = time.time()
            # --- select mating pool ---
            mat_sel = np.argsort(init_score)[0:int(self.mat_size)]
            mat_pop = init_pop[:, mat_sel]
            mat_scores = init_scores[:, mat_sel]
            mat_score = init_score[mat_sel]

            tb = time.time()
            print(f"Select mating population (ms): {1e3*(tb - ta): .3f}")

            # --- crossover ---
            off_pop = np.zeros((len(self.alloc_matrix_poss), self.off_size))
            # take parents in mating pool and produce offsprings
            # midpoint = np.random.randint(0, len(self.alloc_matrix_poss), int(self.pop_size/2))
            midpoint = int(len(self.alloc_matrix_poss) / 2)
            for j in range(int(self.off_size)):
                # take parents in order
                parent1 = mat_pop[:, int(j % self.mat_size)]
                parent2 = mat_pop[:, int((j+1) % self.mat_size)]
                # make offspring by taking first "half" from parent1, second "half" from parent2
                off_pop[:, j] = np.concatenate((parent1[0:midpoint], parent2[midpoint:]))

            tc = time.time()
            print(f"Crossover (ms): {1e3*(tc - tb): .3f}")

            # --- mutation ---
            # select random positioners to allocate to a different target
            nmutation = int(np.ceil(self.pmutation * len(self.alloc_matrix_poss) / 100))
            gene = np.random.randint(0, len(self.alloc_matrix_poss), int(self.off_size) * nmutation)
            for j in range(int(self.mat_size)):
                for k in range(nmutation):
                    # chose new target from the pool of allowed targets for that positioner
                    new_gene = np.random.choice(self.alloc_matrix_poss['reachable'][gene[j * nmutation + k]])
                    off_pop[gene[j * nmutation + k], j] = new_gene

            td = time.time()
            print(f"Mutation (ms): {1e3*(td - tc): .3f}")

            # --- evaluate offspring - time consuming in this form (0.4s for pop_size=20):
            off_scores = self.ga_evaluate(off_pop)
            off_score = np.sum(off_scores, axis=0)

            te = time.time()
            print(f"Evaluate offsprings (ms): {1e3*(te - td): .3f}")

            # --- combine parents and offsprings in new population ---
            init_pop = np.concatenate((mat_pop, off_pop), axis=1)
            init_score = np.concatenate((mat_score, off_score))
            init_scores = np.concatenate((mat_scores, off_scores), axis=1)
            # find best chromosome for plot
            self.solution = init_pop[:, np.argmin(init_score)]
            self.best[i + 1] = np.min(init_score)
            best_scores = init_scores[:, np.argmin(init_score)]
            self.nalloc[i + 1] = len(best_scores[best_scores != 666])

            tf = time.time()
            print(f"Combine offsprings and parent, plot (ms): {1e3 * (tf - te): .3f}")

        t3 = time.time()
        print(f"Generations loops (ms): {1e3 * (t3 - t2)}")
        print(f"Total GA (ms): {1e3 * (t3 - t0)}")

        # Get some info about final solution
        nTar = np.array([len(self.alloc_matrix_poss['reachable'][i]) for i in range(len(self.alloc_matrix_poss))])
        print("Allocatable fibers: ", len(nTar[nTar > 0]))
        print("Allocated fibers: ", self.nalloc[-1])
        print(np.sum(best_scores)/len(nTar[nTar > 0]))

    def ga_evaluate(self, chromosomes):
        """
        Evaluate the score (fitness) of an array of chromosomes
        :param chromosomes: an np.array (nposs, npop)
        :return: scores
        """
        size = np.shape(chromosomes)
        scores = np.zeros_like(chromosomes)
        # TODO: don't always check if gene is already in use by going through the chromosone top to bottom (be random)
        scores = np.array([np.ravel([666 if chromosomes[i, j] in chromosomes[:i, j]
                                     else self.alloc_matrix_poss['distance'][i][chromosomes[i, j] == self.alloc_matrix_poss['reachable'][i]]
                                     for j in range(size[1])])
                           for i in range(size[0])])

        # TODO: apply priority as scaling factor
        # TODO: apply penalty for wide range of magnitude
        # TODO: print range of magnitude for each solution

        # divide by priority of target
        #   score [i] /= self.targets[chromosome[i]]['priority']
        # apply penalty for wide range of magnitudes

        return scores

    def score_apply_priority(self):
        """
        Divide individual distances by the priority of the allocated target
        :return:
        """

        for i in range(len(self.poss)):
            if self.alloc_matrix_poss['score'][i] is not None:
                print(self.targets['priority'][self.alloc_matrix_poss['allocated_1'][i]])
                self.alloc_matrix_poss['score'][i] /= self.targets['priority'][self.alloc_matrix_poss['allocated_1'][i]]

    def score_apply_faint(self):
        """
        Decrease the score of a faint target observed at high tilt and increase that at low tilt
        :return:
        """

        for i in range(len(self.poss)):
            if self.alloc_matrix_poss['score'][i] is not None:
                if self.targets['mag'][i] < 16:
                    return 1
                else:
                    if self.targets['mag'][i] > 24:
                        return - self.alloc_matrix_poss['allocated_dist'][i] / (patrol / 2) + 2
                    else:
                        return (- self.alloc_matrix_poss['allocated_dist'][i] / (patrol / 2) + 2) * (self.targets['mag'][i] - 16) / 8

    def plot(self, ax1, ax2):

        # Mark each PosS with patrol regions
        cc = ['blue', 'green', 'red', 'black']
        for x, y, t, d in zip(self.poss['xpos'], self.poss['ypos'], range(len(self.poss)), self.poss['dither']):
            circle = Circle((x, y), patrol, color=cc[d], alpha=0.5, fill=False)
            ax1.scatter(self.poss['xpos'][self.poss['dither'] == d], self.poss['ypos'][self.poss['dither'] == d],
                        marker='+', c=cc[d])
            ax1.annotate('P-'+str(t)+'-'+str(d), (x, y))
            ax1.add_artist(circle)

        # Mark each Target
        ax1.scatter(self.targets['xpos'], self.targets['ypos'], marker='*', s=10, c='#880088')
        for x, y, t in zip(self.targets['xpos'], self.targets['ypos'], range(len(self.targets))):
            ax1.annotate('T-'+str(t), (x, y))

        # Mark allocated pairs
        for xp, yp, xt, yt in zip(self.poss['xpos'], self.poss['ypos'],
                                  self.targets['xpos'][self.solution], self.targets['ypos'][self.solution]):
            ax1.plot([xp, xt], [yp, yt], 'r+-')

        ax1.set_aspect(aspect=1.)

        # Second plot
        ax2.plot(self.best)

        return

