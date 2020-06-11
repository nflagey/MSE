from mse_alloc import MseFibAlloc

# Without dithering
alloc_nodith = MseFibAlloc(file='cosmo_target_small.csv',
                           doplot=False, meth='one', iternum=1, dither=False, allocfrac=90, spectro='LR')
alloc_nodith.run_optim()
# iternum_nodith = len(alloc_nodith.all_dist)
# all_dist_nodith = [alloc_nodith.all_dist[i].data for i in range(iternum_nodith)]
