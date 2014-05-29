
from __future__ import division
from WindowedModel import WindowedModel  # @UnresolvedImport
from GenotypeModel import GenotypeModel  # @UnresolvedImport
from LDModel import LDModel  # @UnresolvedImport
from IBD.cIBD import cPairIBD
from IBD.intersection import Interval, IntervalTree  # @UnresolvedImport
from collections import Counter 
from IBDAdmixedModel import get_ibdadmixed_inner_models, run_windowed_model

def get_naive_inner_models(ldmodels, K, phased, g=8):
    gpmodels = []
    for anc1 in range(K):
        for anc2 in range(K):
            gpmodels.append(GenotypeModel(ldmodels[anc1],
                                              ldmodels[anc2],
                                              phased=phased,
                                              g=g))
    return gpmodels

def naivemodel(map_file, beagle_model_files, obs_data, ibs_intervals=None, max_snp_num=1000000000, phased=False, g=8, alphas=None, win_size=250, min_score=0):
    K = len(beagle_model_files)
    if alphas is None:
        alphas = [1 / K] * K
        
    ldmodels = []
    for anc in range(K):
        ldmodels.append(LDModel(map_file,
                                beagle_model_files[anc],
                                anc=anc,
                                alpha=alphas[anc],
                                max_snp_num=max_snp_num))
        
    gmodels = get_naive_inner_models(ldmodels, K, phased=phased, g=g)
    local_anc = {}
    breakpoints = {}
    for ind in range(2):
        winmodel = WindowedModel(gmodels, max_snp_num, win_size)
        winmodel.alloc_mem()
        winmodel.calc_ems_probs(obs_data.get_genotype(ind))
        winmodel.viterbi_decoding()
        path = winmodel.get_viterbi_path_models()
        local_anc[ind] = []
        curr_anc1 = -1
        curr_anc2 = -1
        prev_anc1 = -1
        prev_anc2 = -1
        breakpoints[ind] = set()
        win_idx = 0
        for inner_model in path:
            curr_anc1 = min(inner_model._m1._anc, inner_model._m2._anc)
            curr_anc2 = max(inner_model._m1._anc, inner_model._m2._anc)
            if curr_anc1 != prev_anc1 or curr_anc2 != prev_anc2:
                breakpoints[ind].add(win_idx)
            local_anc[ind].append((curr_anc1,curr_anc2))
            prev_anc1 = curr_anc1
            prev_anc2 = curr_anc2
            win_idx += 1
        breakpoints[ind].add(win_idx)
    
    breakpoints_all = breakpoints[0] | breakpoints[1]
    lod_scores = {}
    pairIBD = cPairIBD()
    for interval_idx in range(len(breakpoints_all)-1):
        interval_start = breakpoints_all[interval_idx]
        interval_end = breakpoints_all[interval_idx+1]
        if local_anc[0][interval_start][0] != local_anc[1][interval_start][0] and \
        local_anc[0][interval_start][1] != local_anc[1][interval_start][1] and \
        local_anc[0][interval_start][1] != local_anc[1][interval_start][0] and \
        local_anc[0][interval_start][0] != local_anc[1][interval_start][1]:
            # ancestry in both chromosomes is different. can't be  IBD - just continue
            continue
        else:
            # at least one pair of chromosomes shares the same ancestry. use the most common ancestry for the LD models
            c = Counter(local_anc[0][win_idx][0],
                        local_anc[0][win_idx][1],
                        local_anc[1][win_idx][0],
                        local_anc[1][win_idx][1])
            common_anc = c.most_common(1)[0]
            
            gpmodels_noibd = get_ibdadmixed_inner_models([ldmodels[common_anc]], K=1, ibd=0, phased=phased)
            gpmodels_ibd = get_ibdadmixed_inner_models([ldmodels[common_anc]], K=1, ibd=1, phased=phased)
            
            start_snp = interval_start * win_size
            snp_num = (interval_end - interval_start) * win_size
            interval = (start_snp, start_snp + snp_num) 
            sliced_gpmodels_noibd = [gpmodel.slice_from_model(start_snp,snp_num) for gpmodel in gpmodels_noibd]
            sliced_gpmodels_ibd = [gpmodel.slice_from_model(start_snp,snp_num) for gpmodel in gpmodels_ibd]
            sliced_obs_data = obs_data.get_slice(start_snp,snp_num)
            winmodel_noibd = run_windowed_model(sliced_gpmodels_noibd, sliced_obs_data, sliced_obs_data._snp_num, win_size)
            winmodel_ibd = run_windowed_model(sliced_gpmodels_ibd, sliced_obs_data, sliced_obs_data._snp_num, win_size)
        
            lod_scores[interval] = winmodel_ibd.compare(winmodel_noibd)
            if lod_scores[interval] > min_score:
                pairIBD.add_interval(interval[0],interval[1],lod_scores[interval])         
    
    return (pairIBD,lod_scores)
    
    
    
    
            
