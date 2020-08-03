# from mse_alloc_ga_basic import MseFibAlloc
from mse_alloc_ga import MseFibAlloc

alloc_nodith = MseFibAlloc(file='cosmo_targets_large_new.csv', doplot=False, meth='one', iternum=1, dither=False,
                           allocfrac=90, spectro='LR',
                           init_size=20, mat_size=10, off_size=10, ngeneration=10, pmutation=5)
alloc_nodith.ga_loop()

# Mutation rate:
#   - higher rate means big jumps in fitness with plateaux in between
#   - lower rate means smoother decrease of fitness
#   - final result does not seem to be affected much

# Ngeneration and Off_size are making the code slower (the bottleneck is the evaluation of the offsprings)

# init_size=10, mat_size=10, off_size=30, ngeneration=150, pmutation=1
# --> 2219 out of 2825 allocations (but 215 seconds to do it and fitness was still decreasing steadily)

# init_size=5, mat_size=5, off_size=5, ngeneration=150, pmutation=1
# --> 2135 out of 2825 allocations (in 33 seconds! but fitness was still decreasing steadily)

# init_size=5, mat_size=5, off_size=5, ngeneration=500, pmutation=2
# --> 2261 out of 2825 allocations (in 110 seconds but fitness was still decreasing slowly)

# init_size=5, mat_size=5, off_size=5, ngeneration=1000, pmutation=1
# --> 2312 out of 2825 allocations (in 233 seconds but fitness was still decreasing slowly)

# computing all possible distances was only 100ms while the GAevaluate is typically 200+ ms ...