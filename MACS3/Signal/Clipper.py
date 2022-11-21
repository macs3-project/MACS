import sys
import warnings
import numpy as np
# import pandas as pd
from itertools import combinations
from typing import Hashable, List

# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
# stats = importr('stats')


def match_arg(arg, choices, arg_name):
    if arg in choices:
        return arg
    else:
        sys.exit("arguement \"{}\" should be one of ".format(arg_name) + ", ".join([f"\"{choice}\"" for choice in choices]))


def match(a: List[Hashable], b: List[Hashable]) -> List[int]:
    b_dict = {x: i for i, x in enumerate(b)}
    return [b_dict.get(x, np.nan) for x in a]

def clipper1sided(score_exp, 
                  score_back,
                  regionLength,
                  FDR = 0.05, 
                  ifuseknockoff = None,
                  nknockoff = None,
                  contrastScore_method = None,
                  importanceScore_method = "diff",
                  FDR_control_method = None,
                  ifpowerful = True,
                  seed = 12345):

    score_exp = np.atleast_2d(score_exp)
    score_back = np.atleast_2d(score_back)

    if np.any(score_exp < 0) or np.any(score_back < 0):
        shift = min(np.min(score_exp[~np.isnan(score_exp)]), np.min(score_back[~np.isnan(score_back)]))
        score_exp -= shift
        score_back -= shift
  
    r1 = np.shape(score_exp)[1]
    r2 = np.shape(score_back)[1]

    if np.shape(score_exp)[0] != np.shape(score_back)[0]:
        sys.exit("score_exp and score_back must have the same number of rows (features)")
    
    if ifuseknockoff is None:
        ifuseknockoff = (r1 != r2)

    if not ifuseknockoff:
        if r1 != r2:
            warnings.warn("Caution: no knockoffs are constructed when the numbers of replicates are different; FDR control is not guaranteed")
        if FDR_control_method == "GZ":
            warnings.warn("FDR_control_method cannot be GZ when no knockoffs are constructed. Switching to BC.")
        FDR_control_method = "BC"
        warnings.simplefilter("ignore")
        re = clipper_1sided_woknockoff(score_exp=score_exp, 
                                    score_back=score_back, 
                                    r1=r1, r2=r2, 
                                    FDR=FDR,
                                    importanceScore_method=importanceScore_method,
                                    FDR_control_method=FDR_control_method,
                                    regionLength=regionLength)
    
    if ifuseknockoff:
        if r1 == 1 and r2 == 1:
            sys.exit("Cannot generate knockoffs when both score_exp and score_back have one column. Please rerun clipper1sided by setting ifuseknockoff = F")
        nknockoff_max = len(list(combinations(range(r1+r2), r1)))/2 - 1 if r1 == r2 else len(list(combinations(range(r1+r2), r1))) - 1
        if nknockoff is not None:
            if nknockoff > nknockoff_max or (not isinstance(nknockoff, int)) or nknockoff < 1:
                warnings.warn("nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the exhaustive knockoffs instead.")
            nknockoff = min(nknockoff, nknockoff_max)
        else:
            warnings.warn(f"nknockoff is not supplied; generate the maximum number of knockoffs: {nknockoff_max}")
            nknockoff = nknockoff_max if contrastScore_method == "diff" else 1
        warnings.simplefilter("ignore")
        re = clipper_1sided_wknockoff(score_exp=score_exp,
                                    score_back=score_back,
                                    r1=r1, r2=r2,
                                    FDR=FDR,
                                    importanceScore_method=importanceScore_method,
                                    contrastScore_method= contrastScore_method,
                                    nknockoff=nknockoff,
                                    nknockoff_max=nknockoff_max,
                                    seed=seed)
    
    if ifpowerful and FDR_control_method != "BH":
        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))
        if np.any(FDR_nodisc):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using FDR control method {FDR_control_method}; switching to BH...")
            re_clipperbh = clipper_BH(contrastScore=re["contrastScore"], nknockoff=nknockoff, FDR=FDR[FDR_nodisc])
            re["results"][FDR_nodisc] = re_clipperbh
        
    return re

def clipper2sided(score_exp, 
                  score_back, 
                  FDR = 0.05, 
                  ifuseknockoff = None,
                  nknockoff = None,
                  contrastScore_method = "max",
                  importanceScore_method = "diff",
                  FDR_control_method = "GZ",
                  ifpowerful = True,
                  seed = 12345):
    score_exp = np.atleast_2d(score_exp)
    score_back = np.atleast_2d(score_back)

    r1 = np.shape(score_exp)[1]
    r2 = np.shape(score_back)[1]

    if r1 == 1 and r2 == 1:
        sys.exit("clipper is not yet able to perform two sided identification when either condition has one replicate")
  
    nknockoff_max = (len(list(combinations(range(r1+r2), r1)))/2 - 1) if r1 == r2 else (len(list(combinations(range(r1+r2), r1))) - 1)
    nknockoff_max = int(min(nknockoff_max, 200))
    if nknockoff is not None:
        if nknockoff > nknockoff_max or (not isinstance(nknockoff, int)) or nknockoff < 1:
            warnings.warn("nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the maximal number of knockoffs instead.")
        nknockoff = min(nknockoff, nknockoff_max)
    else:
        nknockoff = min(nknockoff_max, 50) if contrastScore_method == "diff" else 1

    knockoffidx = generate_knockoffidx(r1=r1, r2=r2, 
                                        nknockoff=nknockoff, 
                                        nknockoff_max=nknockoff_max,
                                        seed=seed)

    kappatau_ls = compute_taukappa(score_exp=score_exp,
                                    score_back=score_back,
                                    r1=r1, r2=r2,
                                    if2sided=True,
                                    knockoffidx=knockoffidx,
                                    importanceScore_method=importanceScore_method,
                                    contrastScore_method=contrastScore_method)
    re = clipper_GZ(tau=kappatau_ls["tau"], kappa=kappatau_ls["kappa"],
                    nknockoff=nknockoff, FDR=FDR)
    re = {"knockoffidx": knockoffidx,
            "importanceScore_method": importanceScore_method,
            "importanceScore": kappatau_ls["importanceScore"],
            "contrastScore_method": contrastScore_method,
            "contrastScore": (2*kappatau_ls["kappa"] - 1) * abs(kappatau_ls["tau"]),
            "results": re}

    if ifpowerful and FDR_control_method != "BH":
        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))

        if np.any(FDR_nodisc & (contrastScore_method == "max")):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using FDR control method {FDR_control_method}; switching to BH...")
            re_clipperbh = clipper_BH(contrastScore=re["contrastScore"], nknockoff=nknockoff, FDR=FDR[FDR_nodisc])
            re["results"][FDR_nodisc] = re_clipperbh
    
    return re

def clipper_1sided_woknockoff(score_exp,
                              score_back,
                              r1, r2,
                              importanceScore_method,
                              FDR_control_method,
                              regionLength,
                              FDR = 0.05,
                              aggregation_method = "mean"):
    if r1 > 1:
        score_exp = aggregate_clipper(score=score_exp, aggregation_method=aggregation_method)
    if r2 > 1:
        score_back = aggregate_clipper(score=score_back, aggregation_method=aggregation_method)
    contrastscore = compute_importanceScore_wsinglerep(score_exp=score_exp,
                                                     score_back=score_back,
                                                     importanceScore_method=importanceScore_method)
    if FDR_control_method == "BC":
        warnings.simplefilter("ignore")
        re = clipper_BC(contrastScore=contrastscore, FDR=FDR, regionLength=regionLength)
  
    if FDR_control_method == "BH":
        re = clipper_BH(contrastscore=contrastscore, FDR=FDR)
  
    re = {"importanceScore": contrastscore,
            "importanceScore_method": importanceScore_method,
            "contrastScore": contrastscore,
            "contrastScore_method": importanceScore_method,
            "results": re}
    return re

def clipper_1sided_wknockoff(score_exp,
                             score_back,
                             r1, r2,
                             importanceScore_method,
                             contrastScore_method,
                             nknockoff,
                             nknockoff_max,
                             seed,
                             FDR = 0.05):
    knockoffidx = generate_knockoffidx(r1=r1, r2=r2, 
                                        nknockoff=nknockoff, 
                                        nknockoff_max=nknockoff_max,
                                        seed=seed)
    kappatau_ls = compute_taukappa(score_exp=score_exp,
                                    score_back=score_back,
                                    r1=r1, r2=r2,
                                    if2sided=True,
                                    knockoffidx=knockoffidx,
                                    importanceScore_method=importanceScore_method,
                                    contrastScore_method=contrastScore_method)
    re = clipper_GZ(tau=kappatau_ls["tau"], kappa=kappatau_ls["kappa"], nknockoff=nknockoff, FDR=FDR)
    re = {"knockoffidx": knockoffidx,
            "importanceScore_method": importanceScore_method,
            "importanceScore": kappatau_ls["importanceScore"],
            "contrastScore_method": contrastScore_method,
            "contrastScore": (2*kappatau_ls["kappa"]-1)*abs(kappatau_ls["tau"]),
            "results": re}
    return re


def aggregate_clipper(score, aggregation_method):
    if aggregation_method == "mean":
        score_single = np.array(list(map(np.nanmean, score)))
    if aggregation_method == "median":
        score_single = np.array(list(map(np.nanmedian, score)))
    return score_single

def compute_importanceScore_wsinglerep(score_exp, score_back, importanceScore_method):
    if importanceScore_method == "diff":
        contrastScore = score_exp - score_back
    if importanceScore_method == "max":
        contrastScore = np.maximum(score_exp, score_back)
    return np.array(contrastScore).flatten()

# def clipper_BC(contrastScore, FDR):
#     contrastScore[np.isnan(contrastScore)] = 0
#     c_abs = np.sort(np.unique(np.abs(contrastScore[contrastScore != 0])))
#     i = 0
#     emp_fdp = np.repeat(np.nan, len(c_abs))
#     emp_fdp[0] = 1
#
#     print('contrastScore:')
#     print(contrastScore)
#     # b[np.where(a == value)]
#
#     while (i < len(c_abs)):
#         t = c_abs[i]
#         emp_fdp[i] = min((1 + np.sum(contrastScore <= -t))/np.sum(contrastScore >= t), 1)
#         if i >= 1:
#             emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])
#         i += 1
#
#     c_abs = c_abs[~np.isnan(emp_fdp)]
#     emp_fdp = emp_fdp[~np.isnan(emp_fdp)]
#
#     q_idx = match(contrastScore, c_abs)
#     q = np.array(list(map(lambda k: emp_fdp[k] if k == k else 1, q_idx)))
#
#     def temp(FDR_i):
#         try:
#             thre = c_abs[np.min(np.where(emp_fdp <= FDR_i))]
#         except:
#             thre = np.nan
#         return {"FDR": FDR_i,
#                     "FDR_control": "BC",
#                     "thre": thre,
#                     "q": q,
#                     "discovery": np.where(contrastScore >= thre)[0]
#                     }
#
#     re = np.array(list(map(temp, FDR)))
#     return re

#################################### Add length information (by Saidi) #################################################
def clipper_BC(contrastScore, FDR, regionLength):
    contrastScore[np.isnan(contrastScore)] = 0
    c_abs = np.sort(np.unique(np.abs(contrastScore[contrastScore != 0])))
    i = 0
    emp_fdp = np.repeat(np.nan, len(c_abs))
    emp_fdp[0] = 1

    # print('contrastScore:')
    # print(contrastScore)
    # b[np.where(a == value)]
    # print('regionLength')
    # print(regionLength)
    while (i < len(c_abs)):
        t = c_abs[i]

        # emp_fdp[i] = min((1 + np.sum(contrastScore <= -t)) / np.sum(contrastScore >= t), 1)
        ### Add the length information by Saidi##################
        emp_fdp[i] = min((1 + sum(regionLength[np.where(contrastScore <= -t)])) / sum(regionLength[np.where(contrastScore >= t)]), 1)
        # print(emp_fdp[i])
        if i >= 1:
            emp_fdp[i] = min(emp_fdp[i], emp_fdp[i - 1])
        i += 1
    np.savetxt("/home/saidi/Documents/Postdoc/Project_postdoc/Project1_Clipper_AddOn/MAC3/experiment/output/Clipper_added/syn1/debug_output/emp_fdp",
        np.array(emp_fdp))

    c_abs = c_abs[~np.isnan(emp_fdp)]
    emp_fdp = emp_fdp[~np.isnan(emp_fdp)]

    q_idx = match(contrastScore, c_abs)
    q = np.array(list(map(lambda k: emp_fdp[k] if k == k else 1, q_idx)))

    def temp(FDR_i):
        try:
            thre = c_abs[np.min(np.where(emp_fdp <= FDR_i))]
        except:
            thre = np.nan
        return {"FDR": FDR_i,
                "FDR_control": "BC",
                "thre": thre,
                "q": q,
                "discovery": np.where(contrastScore >= thre)[0]
                }

    re = np.array(list(map(temp, FDR)))
    return re


def clipper_BH(contrastScore, FDR, nknockoff = None):
    if isinstance(contrastScore, dict):
            n = len(np.atleast_1d(contrastScore[1]))
            kappa = contrastScore["kappa"]
            tau = contrastScore["tau"]
            idx_na = np.isnan(tau) | np.isnan(kappa)
            tau = tau[~idx_na]
            kappa = kappa[~idx_na]
            def temp(i):
                t1 = np.array((~kappa) & tau >= tau[i])
                t1 = t1[~np.isnan(t1)]
                t2 = np.array(~kappa)
                t2 = t2[~np.isnan(t2)]
                return sum(t1) / sum(t2) * nknockoff / (nknockoff + 1)
            pval = np.array(list(map(temp, range(n))))
    else:
            contrastScore = np.atleast_1d(contrastScore)
            n = len(contrastScore)
            idx_na = np.isnan(contrastScore)
            contrastScore_nomiss = contrastScore[~idx_na]
            cs_negative = contrastScore_nomiss[contrastScore_nomiss < 0]
            cs_null = np.concatenate((cs_negative, -cs_negative))
            pval = np.array(list(map(lambda x: np.mean(x <= cs_null), contrastScore_nomiss)))

    # qvalue = np.array(stats.p_adjust(FloatVector(pval), method = 'BH'))
    qvalue = pval
    re = np.array(list(map(lambda FDR_i: 
                            {"FDR": FDR,
                                "FDR_control": "BH",
                                "discovery": np.arange(n)[~idx_na][np.where(qvalue <= FDR_i)[0]],
                                "q": qvalue},FDR)))
    return re


def generate_knockoffidx(r1, r2, nknockoff, nknockoff_max, seed):

    np.random.seed(seed)
    if nknockoff_max == 200:
            knockoffidx = [np.nan]*nknockoff
            i_knockoff = 0
            while (i_knockoff < nknockoff):
                temp = np.random.choice(r1+r2, r1, replace=False)
                if np.any(~np.isin(temp, np.arange(r1))) and np.any(np.isin(temp, np.arange(r1))):
                        knockoffidx[i_knockoff] = temp
                        i_knockoff += 1
                else:
                        continue
    else:
        combination_all = np.array([list(x) for x in combinations(range(r1 + r2), r1)], dtype=object)

        if r1 == r2:
            combination_all = combination_all[:len(combination_all)//2][1:]
        else:
            combination_all = combination_all[1:]

        knockoffidx = combination_all[np.random.choice(nknockoff_max, size=nknockoff, replace=False)].tolist()

    return knockoffidx

def compute_taukappa(score_exp, score_back, r1, r2, if2sided,
                     knockoffidx, importanceScore_method, contrastScore_method):

    perm_idx = knockoffidx.copy()
    perm_idx.insert(0, list(range(r1)))

    score_tot = np.concatenate((score_exp, score_back), axis=1)

    def imp_ls_apply(x):
        se = score_tot[:, x]
        sb = score_tot[:, np.setdiff1d(np.arange(r1+r2), x)]
        se = aggregate_clipper(se, aggregation_method='mean')
        sb = aggregate_clipper(sb, aggregation_method='mean')
        imp = compute_importanceScore_wsinglerep(se, sb, importanceScore_method=importanceScore_method)
        if if2sided: imp = np.abs(imp)
        return imp
    imp_ls = np.array(list(map(imp_ls_apply, perm_idx)), dtype=object).T
    def kappatau_ls_apply(x):
        kappa = ~np.any(x[1:] == np.max(x))
        if len(np.atleast_1d(kappa)) == 0:
            kappa = np.nan
        x_sorted = np.sort(x)[::-1]
        if contrastScore_method == "diff":
            tau = x_sorted[0] - x_sorted[1]
        if contrastScore_method == "max":
            tau = x_sorted[0]
        return {"kappa": kappa, "tau": tau}
    kappatau_ls = list(map(kappatau_ls_apply, imp_ls))

    re = {"importanceScore": imp_ls,
            "kappa": np.array([kappatau["kappa"] for kappatau in kappatau_ls]),
            "tau": np.array([kappatau["tau"] for kappatau in kappatau_ls])}
    return re

def clipper_GZ(tau, kappa, nknockoff, FDR):

    contrastScore = (2 * kappa - 1) * np.abs(tau)
    contrastScore[np.isnan(contrastScore)] = 0
    c_abs = np.abs(contrastScore[contrastScore != 0])
    c_abs = np.sort(np.unique(c_abs))

    i = 0
    emp_fdp = np.repeat(np.nan, len(c_abs))
    emp_fdp[0] = 1
    while (i < len(c_abs)):
            t = c_abs[i]
            emp_fdp[i] = min((1/nknockoff + 1/nknockoff * np.sum(contrastScore <= -t))/np.sum(contrastScore >= t), 1)
            if i >= 1: emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])
            i += 1
            
    c_abs = c_abs[~np.isnan(emp_fdp)]
    emp_fdp = emp_fdp[~np.isnan(emp_fdp)]
    q_idx = match(contrastScore, c_abs)
    q = np.array(list(map(lambda k: emp_fdp[k] if k == k else 1, q_idx)))

    def temp(FDR_i):
        try:
            thre = c_abs[np.min(np.where(emp_fdp <= FDR_i))]
        except:
            thre = np.nan
        return {"FDR": FDR_i,
                "FDR_control": "BC",
                "thre": thre,
                "q": q,
                "discovery": np.where(contrastScore >= thre)[0]
                }
    re = np.array(list(map(temp, FDR)))
    return re


def clipper(score_exp, score_back, analysis, regionLength, FDR = 0.05,
            procedure = None, contrast_score = None,
            num_permutation = None, seed = 12345):

    analysis = match_arg(analysis, ["differential", "enrichment"], "analysis")

    score_exp = np.atleast_2d(score_exp)
    score_back = np.atleast_2d(score_back)

    FDR = np.atleast_1d(FDR)
    if (procedure is None):
        procedure = "GZ" if analysis == "differential" else "BC"
    else:
        procedure = match_arg(procedure, ["BC", "aBH", "GZ"], "procedure")

    if (analysis == "differential"):
        contrast_score = match_arg(contrast_score, ["diff", "max"], "contrast_score") if contrast_score is not None else "max"

        re = clipper2sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                        nknockoff=num_permutation, 
                        contrastScore_method= contrast_score, importanceScore_method="diff",
                        FDR_control_method=procedure, ifpowerful=False, seed=seed)
        
        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))
        if np.any(np.array(FDR_nodisc) & (contrast_score == "max")):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold. ")

    elif (analysis == "enrichment"):

        if (np.shape(score_exp)[1] != np.shape(score_back)[1]):
            procedure = "GZ"
        if (contrast_score is None):
            contrast_score = "diff" if procedure == "BC" else "max"
        else:
            contrast_score = match_arg(contrast_score, ["diff", "max"], "contrast_score")
        if procedure == "BC":
            re = clipper1sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                            nknockoff=num_permutation,
                            importanceScore_method=contrast_score,
                            FDR_control_method=procedure, ifpowerful=False, seed=seed,regionLength=regionLength)
        if procedure == "GZ":
            re = clipper1sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                            nknockoff=num_permutation,
                            contrastScore_method=contrast_score,
                            FDR_control_method=procedure, ifpowerful=False, seed=seed,regionLength=regionLength)

        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))
        if np.any(np.array(FDR_nodisc) & (procedure != 'aBH')):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold. ")

    contrast_score_value = re["contrastScore"]
    thre = []
    discoveries = []
    results = re["results"]
    for result in results:
        try:
            thre.append(result["thre"])
            discoveries.append(result["discovery"])
        except:
            continue
    thre = np.array(thre)
    discoveries = np.array(discoveries)
    q = results[0]["q"]
    re = {"contrast.score": contrast_score, 
            "contrast.score.value": contrast_score_value,
            "FDR": FDR,
            "contrast.score.thre": thre,
            "discoveries": discoveries,
            "q": q}

    print('thre:')
    print(thre)

    return re