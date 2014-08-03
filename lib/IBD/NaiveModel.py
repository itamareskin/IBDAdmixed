from __future__ import division
from WindowedModel import WindowedModel  # @UnresolvedImport
from GenotypeModel import GenotypeModel  # @UnresolvedImport
from GenotypePairModel import GenotypePairModel  # @UnresolvedImport
from LDModel import LDModel  # @UnresolvedImport
from IBD.cIBD import cPairIBD
from collections import Counter
from IBDAdmixedModel import run_windowed_model
from collections import OrderedDict

def get_naive_inner_models(ldmodels, K, phased, g=8):
    gpmodels = []
    for anc1 in range(K):
        for anc2 in range(K):
            gpmodels.append(GenotypeModel(m1=ldmodels[anc1],
                                          m2=ldmodels[anc2],
                                          phased=phased,
                                          g=g))
    return gpmodels

def get_naive_inner_models_second_stage(ldmodel, phased, ibd, g):
    return GenotypePairModel(ldmodel,
                              ldmodel,
                              ldmodel,
                              ldmodel,
                              phased=phased,
                              ibd=ibd,
                              g=g)

def get_ancestry_breakpoints(winmodel, obs_data):
    local_anc = {}
    breakpoints = {}
    print "Running local ancestry inference..."
    for ind in range(2):
        print "Running local ancestry inference of individual: " + str(ind)
        print "running viterbi decoding..."
        winmodel.viterbi_decoding(obs_data.get_genotype(ind))
        print "running viterbi backtracking..."
        path = winmodel.get_viterbi_path_models()
        local_anc[ind] = []
        prev_anc1 = -1
        prev_anc2 = -1
        breakpoints[ind] = set()
        win_idx = 0
        print "Finding breakpoints..."
        for inner_model in path:
            curr_anc1 = min(inner_model._m1._anc, inner_model._m2._anc)
            curr_anc2 = max(inner_model._m1._anc, inner_model._m2._anc)
            if curr_anc1 != prev_anc1 or curr_anc2 != prev_anc2:
                breakpoints[ind].add(win_idx)
            local_anc[ind].append((curr_anc1, curr_anc2))
            prev_anc1 = curr_anc1
            prev_anc2 = curr_anc2
            win_idx += 1
        breakpoints[ind].add(win_idx)
    print "Finished finiding all breakpoints."
    return (sorted(list(breakpoints[0] | breakpoints[1])), local_anc)

def naivemodel(map_file, beagle_model_files, obs_data, ibs_intervals=None, max_snp_num=1000000000, phased=False, g=8,
               alphas=None, win_size=250, min_score=0, offset=25):
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
    snp_num = ldmodels[0]._snp_num
    gmodels = get_naive_inner_models(ldmodels, K, phased=phased, g=g)
    winmodel = WindowedModel(gmodels, snp_num, win_size)
    (breakpoints, local_anc) = get_ancestry_breakpoints(winmodel, obs_data)

    pairIBD = cPairIBD()
    lod_scores = {}
    print "Running Naive model between all ancestry breakpoints..."
    for interval_idx in range(len(breakpoints) - 1):
        interval_start = breakpoints[interval_idx]
        interval_end = breakpoints[interval_idx + 1]
        interval = (interval_start * win_size, interval_end * win_size)
        print "working on interval: " + str(interval)
        max_offset = win_size
        lod_scores[interval] = {}
        for offset_idx in range(int(max_offset/offset)):
            print "working on offset: " + str(offset_idx)
            start_snp = interval[0] + offset_idx*offset
            snp_num = interval[1] - interval[0]
            if local_anc[0][interval_start][0] == local_anc[1][interval_start][0] or \
                            local_anc[0][interval_start][1] == local_anc[1][interval_start][1] or \
                            local_anc[0][interval_start][1] == local_anc[1][interval_start][0] or \
                            local_anc[0][interval_start][0] == local_anc[1][interval_start][1]:
                # at least one pair of chromosomes shares the same ancestry. use the most common ancestry for the LD models
                print "Running IBD model..."
                c = Counter([local_anc[0][interval_start][0],
                             local_anc[0][interval_start][1],
                             local_anc[1][interval_start][0],
                             local_anc[1][interval_start][1]])
                common_anc = c.most_common(1)[0][0]

                gpmodel_noibd = get_naive_inner_models_second_stage(ldmodels[common_anc],phased,0,g)
                gpmodel_ibd = get_naive_inner_models_second_stage(ldmodels[common_anc],phased,1,g)

                sliced_gpmodel_noibd = gpmodel_noibd.slice_from_model(start_snp, snp_num)
                sliced_gpmodel_ibd = gpmodel_ibd.slice_from_model(start_snp, snp_num)
                sliced_obs_data = obs_data.get_slice(start_snp, snp_num)
                winmodel = run_windowed_model([[sliced_gpmodel_noibd], [sliced_gpmodel_ibd]], sliced_obs_data, sliced_obs_data._snp_num, win_size)
                win_scores = winmodel.compare_partitions(1,0)

                for win in win_scores.keys():
                    lod_scores[interval][start_snp+win[0],start_snp+win[1]] = win_scores[win]
                    if win_scores[win] > min_score:
                        pairIBD.add_interval(start_snp+win[0],start_snp+win[1],win_scores[win])
        lod_scores[interval] = OrderedDict(sorted(lod_scores[interval].items(), key=lambda t: t[0][0]))
    lod_scores = OrderedDict(sorted(lod_scores.items(), key=lambda t: t[0][0]))
    if len(pairIBD.to_list()) > 0:
        pairIBD.merge_intervals(max_val=True,merge_diff_vals=True)
    return (pairIBD, lod_scores)
    
    
    
    
            
