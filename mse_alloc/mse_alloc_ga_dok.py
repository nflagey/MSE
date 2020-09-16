#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, vstack
from astropy import wcs
import time
from scipy.spatial import cKDTree
import random
import matplotlib.pyplot as plt
from numpy.random import default_rng

# Random number generator
rng = default_rng()

# Global variables
plate_scale = 106.7e-3  # microns per arcsec
patrol = 1.24 * 7.77 / plate_scale  # in arcsec
bigR = 7.77 / plate_scale  # circum-radius
smallR = np.sqrt(3) / 2 * bigR  # in-radius

# About fitness
score_none = 1.e6  # distance for unallocated fiber's score
fitness_method = 'distance'  # 'alloc' for number of allocation, 'distance' for distance

# GA parameters
pop_size = 10  # size of population
elite_pct = 20.  # percent of best parents to keep for next generation
mating_pct = 50.  # percent of best parents allowed to mate
mutation_pct = .5  # percent of genes that mutate after crossover
max_generation = 100  # maximum number of new generations
max_time = 600  # maximum time for the whole process


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Trying a different approach for the class/methods/functions
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

class FiberToTargetAllocation(object):
    def __init__(self, ind_tgt_near_poss, n_targets, distances_dok, targets, first_gen=False, distances_coo=None):
        self.ind_tgt_near_poss = ind_tgt_near_poss
        self.n_targets = n_targets
        self.distances_dok = distances_dok
        self.distances_coo = distances_coo
        self.targets = targets
        self.first_gen = first_gen

        # Execute on init:
        self._chromosome = None
        self._fitness = None
        self._fitness_clean = None
        self._distances = None
        self._priorities = None

    @property
    def chromosome(self):
        if getattr(self, '_chromosome', None) is None:
            chromosome = np.zeros(len(self.ind_tgt_near_poss), dtype=np.int32) - 1  # faster than list comprehension

            # TODO: try to reproduce basic algorithm here as the firs generation?
            # pick up a random order
            # go through the fibers
            # pick the closest target (or shortest distance/priority)
            # "nullify" that row and column from distance/ind_tgt_near_poss matrices
            # continue until all possible are done
            # order = np.random.choice(replace=False)



            if self.first_gen:
                # allocate targets to their nearest fibers
                # chromosome = np.argmin(self.distances_arr, axis=0)  # very slow with dense array
                for i, n in enumerate(self.n_targets):
                    if n > 0:
                        column = self.distances_coo.col == i
                        chromosome[i] = self.distances_coo.row[column][np.argmin(self.distances_coo.data[column])]
            else:
                # allocate targets that are the only ones a fiber can reach
                chromosome[self.n_targets == 1] = self.ind_tgt_near_poss[self.n_targets == 1][0]
                # chose randomly for the fibers that can reach more than one target
                a = list(map(lambda x: rng.choice(x), self.ind_tgt_near_poss[self.n_targets > 1]))
                chromosome[self.n_targets > 1] = a

            # Take care of duplicates: remove all but one fiber in those cases
            uniq, counts = np.unique(chromosome, return_counts=True)
            for u, c in zip(uniq, counts):
                # is that target allocated more than once?
                if c > 1 and u != -1:
                    # find all fibers allocated to that target
                    ind_dup = np.ravel(np.argwhere(chromosome == u))
                    # put all to -1
                    chromosome[ind_dup] = -1
                    # randomly chose one to be back to u
                    chromosome[rng.choice(ind_dup)] = u

                    # all solutions below are slower:
                    # pick all but one with random.sample: 1400ms
                    # chromosome[random.sample(list(ind_dup), c-1)] = np.nan
                    # pick all but one with np.random.choice: 2000ms
                    # chromosome[np.random.choice(ind_dup, size=c-1, replace=False)] = np.nan
                    # np.random.shuffle: 1200ms
                    # np.random.shuffle(ind_dup)
                    # chromosome[ind_dup[:-1]] = np.nan

            self._chromosome = chromosome
        return self._chromosome

    @chromosome.setter
    def chromosome(self, value):
        self._chromosome = value

    def mate(self, parent2):
        # Check all probabilities at once (0: mutate, 1: parent1, 2: parent2)
        prob = rng.choice([0, 1, 2], p=[mutation_pct/100, (100-mutation_pct)/200, (100-mutation_pct)/200],
                                size=len(self.ind_tgt_near_poss))
        # Create a third chromosome all at once to store the mutated genes
        mutated = FiberToTargetAllocation(self.ind_tgt_near_poss, self.n_targets, self.distances_dok, self.targets)
        # Now generate offspring
        offspring = FiberToTargetAllocation(self.ind_tgt_near_poss, self.n_targets, self.distances_dok, self.targets)
        offspring.chromosome = mutated.chromosome
        offspring.chromosome[prob == 1] = self.chromosome[prob == 1]
        offspring.chromosome[prob == 2] = parent2.chromosome[prob == 2]

        return offspring

    @property
    def fitness(self):
        if getattr(self, '_fitness', None) is None:
            self._fitness = self.compute_fitness()

        return self._fitness

    @property
    def fitness_clean(self):
        if getattr(self, '_fitness_clean', None) is None:
            self._fitness_clean = self.compute_fitness(score_none=0)

        return self._fitness_clean

    @property
    def priorities(self):
        if getattr(self, '_priorities', None) is None:
            priorities = np.zeros(len(self.chromosome))
            for i, c in enumerate(self.chromosome):  # i is the positioner index, c is the target index
                if c != -1:
                    priorities[i] = (self.targets[c]['priority'] * self.targets[c]['surveypriority'])
            self._priorities = priorities
        return self._priorities

    @property
    def distances(self):
        if getattr(self, '_distances', None) is None:
            distances = np.zeros(len(self.chromosome))
            for i, c in enumerate(self.chromosome):  # i is the positioner index, c is the target index
                if c != -1:
                    distances[i] = self.distances_dok[(c, i)]
            self._distances = distances
        return self._distances

    def compute_fitness(self, score_none=score_none):
        fitness = np.zeros(len(self.chromosome))
        for i, c in enumerate(self.chromosome):  # i is the positioner index, c is the target index
            if c == -1:
                fitness[i] = score_none
            else:
                # DOK matrix is about twice as fast as the record array
                # get distance
                fitness[i] = self.distances_dok[(c, i)]
                # apply priority score
                fitness[i] /= (self.targets[c]['priority'] * self.targets[c]['surveypriority'])

        # Account for range of magnitudes
        # Get magnitudes of allocated targets
        sel = self.chromosome[self.chromosome >= 0]
        umags, gmags, rmags = self.targets['umag'][sel], self.targets['gmag'][sel], self.targets['rmag'][sel]
        imags, zmags = self.targets['imag'][sel], self.targets['zmag'][sel]
        jmags, hmags = self.targets['Jmag'][sel], self.targets['Hmag'][sel]
        # Use standard deviations instead of max-min
        umag_rms, gmag_rms, rmag_rms = np.nanstd(umags), np.nanstd(gmags), np.nanstd(rmags)
        imag_rms, zmag_rms = np.nanstd(imags), np.nanstd(zmags)
        jmag_rms, hmag_rms = np.nanstd(jmags), np.nanstd(hmags)
        # Compute scaling
        scale = (umag_rms + gmag_rms + rmag_rms + imag_rms + zmag_rms + jmag_rms + hmag_rms) / 7. + 1
        # make sure we always maximize
        fitness = -np.sum(fitness) * scale  # we thus want to get as close to 0 as possible here

        return fitness


def create_targets(file):
    # Read targets file
    targets = Table.read('TARGETS/' + file, format='csv')
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
    # Add coordinates (center on mean RA and mean DEC)
    hdr['CRVAL1'] = np.mean(targets['RAJ2000'])
    hdr['CRVAL2'] = np.mean(targets['DECJ2000'])
    # Get pixel coordinates for all targets
    wcs0 = wcs.WCS(hdul[0])
    xy = [wcs0.wcs_world2pix([[targets['RAJ2000'][i], targets['DECJ2000'][i]]], 1)
          for i in range(len(targets['RAJ2000']))]
    # Make those columns
    xpos = Column([xy[i][0][0] for i in range(len(xy))], name='xpos', unit='arcsec')
    ypos = Column([xy[i][0][1] for i in range(len(xy))], name='ypos', unit='arcsec')
    nobs = Column(np.zeros(len(xy)), name='nobs')
    # Add columns to table
    targets.add_columns([xpos, ypos, nobs])
    # Print info about targets
    # print(np.nanmin(targets['umag']), np.nanmax(targets['umag']))
    # print(np.nanmin(targets['gmag']), np.nanmax(targets['gmag']))
    # print(np.nanmin(targets['rmag']), np.nanmax(targets['rmag']))
    # print(np.nanmin(targets['imag']), np.nanmax(targets['imag']))
    # print(np.nanmin(targets['zmag']), np.nanmax(targets['zmag']))
    # print(np.nanmin(targets['Jmag']), np.nanmax(targets['Jmag']))
    # print(np.nanmin(targets['Hmag']), np.nanmax(targets['Hmag']))
    # Store
    return targets


def create_poss(dither=False, spectro='LR'):
    # Read AAO file (coordinates in mm)
    poss = Table.read('POSS/' + 'sphinx.tsv', format='ascii')
    # only select LR or HR
    poss = poss[poss['col3'] == spectro]
    # convert into arcsec
    poss['col1'] /= plate_scale
    poss['col2'] /= plate_scale
    # rename columns
    poss['col1'].name = 'xpos'
    poss['col2'].name = 'ypos'
    poss['col3'].name = 'type'
    # add dither column
    poss['dither'] = 0
    # Dithering?
    if dither:
        print("<h3> You are using Dithering<h3>")
        poss = vstack([poss,
                       Table([poss['xpos'] + bigR / 4, poss['ypos'] + smallR / 2, poss['type'], poss['dither'] + 1],
                             names=('xpos', 'ypos', 'type', 'dither')),
                       Table([poss['xpos'] + 3 * bigR / 4, poss['ypos'] + smallR / 2, poss['type'], poss['dither'] + 2],
                             names=('xpos', 'ypos', 'type', 'dither')),
                       Table([poss['xpos'] + bigR / 2, poss['ypos'], poss['type'], poss['dither'] + 3],
                             names=('xpos', 'ypos', 'type', 'dither'))])

    return poss


def compute_distances(targets, poss):
    # extract XY columns only for each of the LR and HR
    poss_xy = np.array([poss['xpos'], poss['ypos']]).transpose()
    targets_xy = np.array([targets['xpos'], targets['ypos']]).transpose()
    # get kd-trees
    kdt_target = cKDTree(targets_xy)  # fibers to targets
    kdt_poss = cKDTree(poss_xy)  # targets to fibers
    # get nearby targets for each fiber --> this returns a list of lists (super fast!)
    ind_target_near_poss = np.asarray(kdt_poss.query_ball_tree(kdt_target, patrol, p=2))
    n_target = np.array([len(ind_target_near_poss[i]) for i in range(len(poss))])
    # remove fibers that cannot reach any target
    # ind_target_near_poss = ind_target_near_poss[n_target > 0] # Using -1 in the chromosome instead

    # get nearby fibers for each target
    ind_poss_near_target = np.asarray(kdt_target.query_ball_tree(kdt_poss, patrol, p=2))
    n_poss = np.array([len(ind_poss_near_target[i]) for i in range(len(targets))])
    # remove targets that cannot be reached by any fiber
    # ind_poss_near_target = ind_poss_near_target[n_poss > 0] # Using -1 in the chromosome instead

    # Compute distances and put everything in sparse DoK matrix
    dist_target_dok = kdt_target.sparse_distance_matrix(kdt_poss, patrol, p=2.)

    #   access dok-matrix this way: > dist_target[(itarget, jposs)])

    # TODO: what if the distance between a fiber and a target is exactly 0? Will the sparse matrix remove it?
    # TODO: better to keep all targets/fibers and pass them to the GA or only pass those that need to be optimized?

    return dist_target_dok, ind_target_near_poss, ind_poss_near_target, n_target, n_poss


def main(file, t_start):
    # read input file and create targets
    targets = create_targets(file)
    # create PosS
    poss = create_poss(dither=False, spectro='LR')
    # compute distances, figure out reachable targets and available fibers
    distances_dok, ind_tgt_near_poss, ind_poss_near_tgt, n_target, n_poss = compute_distances(targets, poss)
    # distances_arr = distances_dok.toarray()
    distances_coo = distances_dok.tocoo()
    # distances_arr[distances_arr == 0] = score_none

    # Look at some statistics
    #  how many fibers can be allocated to at least one target?
    max_alloc = np.count_nonzero(n_target > 0)
    #  how many fibers can be allocated to exactly one target?
    # print(np.count_nonzero(n_target == 1))
    #  how many possible different chromosome total?
    doubles = (n_target[n_target > 1])
    # print(len(doubles))

    # create initial population and evaluate them automatically
    # t0 = time.time()
    population = []
    for i in range(pop_size):
        x = FiberToTargetAllocation(ind_tgt_near_poss, n_target, distances_dok, targets,
                                    first_gen=True, distances_coo=distances_coo)
        population.append(x)
    # at this stage the chromosomes have not been generated, we only created instances of the class

    # t1 = time.time()
    # print(f"Initial population (ms): {1e3 * (t1 - t0):.3f}")

    # go through generations
    all_n_alloc = np.zeros(max_generation)
    optimized = False
    generation = 0
    while not optimized:
        # t2 = time.time()
        # sort the population in increasing order of fitness score
        population = sorted(population, key=lambda x: x.fitness, reverse=True)
        # to get the fitness, we need the chromosome to be generated so this is when it's happening

        # t3 = time.time()
        # print(f"Evaluate population (ms): {1e3*(t3 - t2): .3f}")

        # Figure out how many fiber/target are paired
        n_alloc = np.count_nonzero(population[0].chromosome >= 0)
        all_n_alloc[generation] = n_alloc
        # Plot
        plt.figure(1)
        plt.subplot(221)
        plt.plot(generation, n_alloc, 'r+')
        plt.subplot(222)
        plt.plot(generation, population[0].fitness, 'bx')
        plt.subplot(223)
        plt.plot(generation, np.mean(population[0].distances[population[0].distances > 0]), 'g*')
        plt.subplot(224)
        plt.plot(generation, np.mean(population[0].priorities[population[0].priorities > 0]), 'bx')

        # Break now if we reached the max number of generations
        generation += 1
        if generation >= max_generation:
            print("Max generation reached")
            optimized = True
            break
        # Break now if we reached the max number of allocation
        if n_alloc == max_alloc:
            print("Max number of allocation reached")
            optimized = True
            break
        # Break now if we reached the max amount of time
        if time.time() - t_start > max_time:
            print("Max optimization time reached")
            optimized = True
            break
        # Break now if we reached a plateau
        if generation > 10:
            if np.min(all_n_alloc[generation-11:generation-1]) > all_n_alloc[generation-1]:
                print("Plateau reached")
                optimized = True
                break

        # prepare for next generation
        next_generation = []

        # elitism: keep elite_pct% of the best parents
        s = int(elite_pct * pop_size / 100)
        next_generation.extend(population[:s])

        # t4 = time.time()
        # print(f"Practicing elitism (ms): {1e3 * (t4 - t3): .3f}")

        # mating: use mating_pct% of the best parents to generate the rest of the next generation
        s = int((100 - elite_pct) * pop_size / 100)
        m = int(mating_pct * pop_size / 100)
        for _ in range(s):
            parent1 = rng.choice(population[:m])
            parent2 = rng.choice(population[:m])
            offspring = parent1.mate(parent2)
            next_generation.append(offspring)

        # t5 = time.time()
        # print(f"Crossover & mutation (ms): {1e3 * (t5 - t4): .3f}")

        population = next_generation

    plt.subplot(221)
    plt.title(f"Number of allocations")
    plt.xlabel("Number of generations")
    plt.subplot(222)
    plt.title(f"Total fitness")
    plt.xlabel("Number of generations")
    plt.subplot(223)
    plt.title(f"Average distance per target")
    plt.xlabel("Number of generations")
    plt.subplot(224)
    plt.title(f"Average priority per target")
    plt.xlabel("Number of generations")
    # compute typical priority score for targets that can be reached
    ref_priorities = targets[n_poss > 0]['priority'] * targets[n_poss > 0]['surveypriority']
    plt.plot([0, generation], [np.mean(ref_priorities), np.mean(ref_priorities)], 'r')
    plt.show()

    return n_alloc


if __name__ == '__main__':
    t_start = time.time()
    file = 'cosmo_targets_large_new.csv'
    n_alloc = main(file, t_start)
    print(f"Final number of allocated pairs: {n_alloc}")
    print(f"Total time (s): {(time.time() - t_start): .3f}")
