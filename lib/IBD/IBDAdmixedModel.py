
from __future__ import division
from WindowedModel import WindowedModel  # @UnresolvedImport
from GenotypePairModel import GenotypePairModel  # @UnresolvedImport
from LDModel import LDModel  # @UnresolvedImport
from IBD.cIBD import cPairIBD
from collections import OrderedDict

def get_ibdadmixed_inner_models(ldmodels, K, phased, ibd, g):
    gpmodels = []
    for anc1 in range(K):
        for anc2 in range(K):
            for anc3 in range(K):
                for anc4 in range(K):
                    gpmodels.append(GenotypePairModel(ldmodels[anc1],
                                                      ldmodels[anc2],
                                                      ldmodels[anc3],
                                                      ldmodels[anc4],
                                                      phased=phased,
                                                      ibd=ibd,
                                                      g=g))
    return gpmodels

def run_windowed_model(gpmodels, obs_data, max_snp_num, win_size):
    winmodel = WindowedModel(gpmodels, max_snp_num, win_size)
    winmodel.posterior_decoding(obs_data, False)
    return winmodel

def ibdadmixed(map_file, beagle_model_files, obs_data, ibs_intervals=None, max_snp_num=1000000000, phased=False, g=8, alphas=None, ibd_trans=None, win_size=250, min_score=0, offset=25):
    K = len(beagle_model_files)
    if alphas is None:
        alphas = [1 / K] * K
    if ibd_trans is None:
        ibd_trans = [1e-5, 1] * K
        
    ldmodels = []
    for anc in range(K):
        ldmodels.append(LDModel(map_file,
                                beagle_model_files[anc],
                                anc=anc,
                                alpha=alphas[anc],
                                t_0_1=ibd_trans[anc * 2],
                                t_1_0=ibd_trans[anc * 2 + 1],
                                max_snp_num=max_snp_num))
        
    gpmodels_noibd = get_ibdadmixed_inner_models(ldmodels, K, phased=phased, ibd=False, g=g)
    gpmodels_ibd = get_ibdadmixed_inner_models(ldmodels, K, phased=phased, ibd=True, g=g)

    pairIBD = cPairIBD()
    lod_scores = {}
    if ibs_intervals is None:
        ibs_intervals = [(0, obs_data._snp_num)]
    print "Running IBDAdmixed on IBS intervals..."
    for interval in ibs_intervals:
        print "working on interval: " + str(interval)
        max_offset = win_size
        lod_scores[interval] = {}
        for offset_idx in range(int(max_offset/offset)):
            print "working on offset: " + str(offset_idx)
            start_snp = interval[0] + offset_idx*offset
            snp_num = interval[1] - interval[0]
            sliced_gpmodels_noibd = [gpmodel.slice_from_model(start_snp,snp_num) for gpmodel in gpmodels_noibd]
            sliced_gpmodels_ibd = [gpmodel.slice_from_model(start_snp,snp_num) for gpmodel in gpmodels_ibd]
            sliced_obs_data = obs_data.get_slice(start_snp,snp_num)
            winmodel = run_windowed_model([sliced_gpmodels_noibd, sliced_gpmodels_ibd], sliced_obs_data, sliced_obs_data._snp_num, win_size)
            win_scores = winmodel.compare_partitions(1,0)
            for win in win_scores.keys():
                lod_scores[interval][start_snp+win[0],start_snp+win[1]] = win_scores[win]
                if win_scores[win] > min_score:
                    pairIBD.add_interval(start_snp+win[0],start_snp+win[1],win_scores[win])
        lod_scores[interval] = OrderedDict(sorted(lod_scores[interval].items(), key=lambda t: t[0][0]))
    lod_scores = OrderedDict(sorted(lod_scores.items(), key=lambda t: t[0][0]))
    pairIBD.merge_intervals(max_val=True,merge_diff_vals=True)
    return (pairIBD, lod_scores)
    
    
    
    
            
