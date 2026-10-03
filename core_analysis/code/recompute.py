"""Independent arithmetic review, using deposited/processed data supplied by authors.

This is a post hoc review, not a replacement for a prospectively frozen protocol.
Requires Python 3.10+, numpy and scipy. See README.md for scope and invocation.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import statistics
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize_scalar
from scipy.stats import t


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_zip(path, name):
    with zipfile.ZipFile(path) as archive:
        return list(csv.DictReader(io.StringIO(archive.read(name).decode("utf-8-sig"))))


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metrics(cases, controls):
    x, y = np.asarray(cases, dtype=float), np.asarray(controls, dtype=float)
    assert len(x) >= 2 and len(y) >= 2
    assert np.isfinite(x).all() and np.isfinite(y).all()
    pair = (x[:, None] > y).astype(float) + 0.5 * (x[:, None] == y)
    auc = float(pair.mean())
    av = pair.mean(axis=1).var(ddof=1) / len(x) + pair.mean(axis=0).var(ddof=1) / len(y)
    ase = math.sqrt(float(av))
    df = len(x) + len(y) - 2
    pooled_var = ((len(x)-1)*x.var(ddof=1) + (len(y)-1)*y.var(ddof=1)) / df
    if pooled_var <= 0:
        raise ValueError("Zero within-group variance; no artificial precision is assigned.")
    d = (x.mean()-y.mean()) / math.sqrt(pooled_var)
    correction = 1 - 3 / (4*(len(x)+len(y))-9)
    g = float(correction*d)
    se = math.sqrt(correction**2 * ((len(x)+len(y))/(len(x)*len(y)) + d*d/(2*df)))
    return dict(N_Psoriasis=len(x), N_Healthy=len(y), AUC=auc, AUC_SE=ase,
                AUC_Lower=max(0, auc-1.96*ase), AUC_Upper=min(1, auc+1.96*ase),
                AUC_IntervalStatus="DEGENERATE_ZERO_SE" if ase == 0 else "NORMAL_APPROXIMATION",
                HedgesG=g, HedgesG_SE=se, HedgesG_Lower=g-1.96*se, HedgesG_Upper=g+1.96*se,
                OrientedPsoriasisMean=float(x.mean()), OrientedHealthyMean=float(y.mean()),
                OrientedMeanDifference=float(x.mean()-y.mean()))


def classification(lo, hi):
    return "POSITIVE" if lo > 0 else "NEGATIVE" if hi < 0 else "CROSSES_ZERO"


def corrected_meta(rows):
    """REML objective minimized over a domain bounded by the effect range.

    Beyond the positive root of (k-1)*tau2^2 =
    k*range(y)^2*(max(v)+tau2), the derivative is nonnegative.
    We use twice that bound (with a floor), inspect a dense grid, and refine every local
    minimum as well as the zero boundary. This avoids the original
    premature upper-bound truncation. These are numerical checks, not a
    new statistical estimator.
    """
    y = np.array([float(row["HedgesG_Aligned"]) for row in rows])
    v = np.array([float(row["HedgesG_SE"])**2 for row in rows])
    k = len(y)
    assert k >= 2 and (v > 0).all()
    def objective(tau):
        w = 1/(v+tau)
        mu = (w*y).sum()/w.sum()
        return float(.5*(np.log(v+tau).sum()+np.log(w.sum())+(w*(y-mu)**2).sum()))
    # A sufficient conservative bound follows from
    # (k-1)/(vmax+tau) >= k*range(y)^2/tau^2.
    a = k*float(np.ptp(y))**2/(k-1)
    bound = max(1.0, (a+math.sqrt(a*a+4*a*float(v.max())))/2)*2
    scale = max(float(v.min())*1e-6, 1e-12)
    grid = np.r_[0., np.geomspace(scale, bound, 801)]
    vals = np.array([objective(q) for q in grid])
    candidates = [(objective(0), 0.)]
    for i in range(1, len(grid)-1):
        if vals[i] <= vals[i-1] and vals[i] <= vals[i+1]:
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]),
                                  method="bounded", options={"xatol":1e-12})
            candidates.append((float(fit.fun), float(fit.x)))
    fit = minimize_scalar(objective, bounds=(0, bound), method="bounded", options={"xatol":1e-12})
    candidates.append((float(fit.fun), float(fit.x)))
    best, tau = min(candidates)
    if objective(0) <= best+1e-10:
        tau = 0.
    w = 1/(v+tau)
    mu = float((w*y).sum()/w.sum())
    se = math.sqrt(1/w.sum())
    qstar = max(1, float((w*(y-mu)**2).sum())/(k-1))
    mhkse = se*math.sqrt(qstar)
    crit = float(t.ppf(.975,k-1))
    wf = 1/v
    q = float((wf*(y-(wf*y).sum()/wf.sum())**2).sum())
    return dict(K=k, Tau2=tau, HedgesG=mu, SE=se, CI_Lower=mu-1.96*se, CI_Upper=mu+1.96*se,
                I2=max(0,(q-k+1)/q*100) if q else 0,
                mHK_Lower=mu-crit*mhkse, mHK_Upper=mu+crit*mhkse,
                NormalPI_Lower=mu-1.96*math.sqrt(tau+se*se), NormalPI_Upper=mu+1.96*math.sqrt(tau+se*se),
                mHK_t_PI_Lower=mu-crit*math.sqrt(tau+mhkse*mhkse),
                mHK_t_PI_Upper=mu+crit*math.sqrt(tau+mhkse*mhkse))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrices", required=True, type=Path)
    parser.add_argument("--stage4a", required=True, type=Path)
    parser.add_argument("--stage4f", required=True, type=Path)
    parser.add_argument("--stage5b", required=True, type=Path)
    parser.add_argument("--holdout", required=True, type=Path)
    parser.add_argument("--exclusions", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    old = read_zip(args.stage4a,"single_gene_per_cohort_metrics_STAGE4A.csv")
    pooled = read_zip(args.stage4f,"Table1_Primary_Marker_Results_STAGE4F.csv")
    holdold = read_zip(args.stage5b,"GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv")
    excludes = None
    if args.exclusions:
        excludes = {r["GSM"].strip() for r in read_csv(args.exclusions)
                    if r.get("SampleExclusionSetID","GSE54456_SHARED_42_OFFICIAL") == "GSE54456_SHARED_42_OFFICIAL"}
        assert len(excludes) == 42
    matrixfiles = sorted(args.matrices.glob("*_primary_gene_matrix_STAGE3C.csv"))
    matrices = {p.name.split("_")[0]:read_csv(p) for p in matrixfiles}
    assert len(matrices) == 7
    sample_audit = []
    for gse, rs in matrices.items():
        assert len({r["GSM"] for r in rs}) == len(rs)
        sample_audit.append(dict(GSE=gse, Samples=len(rs), Genes=len(rs[0])-6,
                                 Psoriasis=sum(r["FinalPhenotype"].startswith("PSORIASIS_LESIONAL_SKIN") for r in rs),
                                 Healthy=sum(r["FinalPhenotype"]=="HEALTHY_CONTROL_SKIN" for r in rs)))
    checks, pending, recomputed = [], [], {}
    for r in old:
        if int(r["ExcludedN"]) and excludes is None:
            pending.append(dict(TestID=r["TestID"], PMID=r["PMID"], Gene=r["GeneSymbol"], GSE=r["GSE"],
                                Reason="Missing official 42-sample exclusion list"))
            continue
        data = [x for x in matrices[r["GSE"]] if not int(r["ExcludedN"]) or x["GSM"] not in excludes]
        ori = r["FrozenOrientation"]
        assert ori in ("UP_IN_PSORIASIS","DOWN_IN_PSORIASIS")
        sign = 1 if ori == "UP_IN_PSORIASIS" else -1
        cases = [sign*float(x[r["GeneSymbol"]]) for x in data if x["FinalPhenotype"].startswith("PSORIASIS_LESIONAL_SKIN")]
        controls = [sign*float(x[r["GeneSymbol"]]) for x in data if x["FinalPhenotype"]=="HEALTHY_CONTROL_SKIN"]
        assert (len(cases),len(controls)) == (int(r["N_Psoriasis"]),int(r["N_Healthy"]))
        m = metrics(cases,controls)
        delta = max(abs(m["AUC"]-float(r["AUC_Directional"])),abs(m["AUC_SE"]-float(r["AUC_SE_DeLong"])),
                    abs(m["HedgesG"]-float(r["HedgesG_Aligned"])),abs(m["HedgesG_SE"]-float(r["HedgesG_SE"])))
        checks.append(dict(TestID=r["TestID"], PMID=r["PMID"], Gene=r["GeneSymbol"], GSE=r["GSE"],
                           **m, MaxAbsoluteDifference=delta, MatchesWithin1e_6=delta<1e-6))
        recomputed[r["TestID"]] = dict(r, HedgesG_Aligned=m["HedgesG"],HedgesG_SE=m["HedgesG_SE"])
    write_csv(args.out/"development_recomputed.csv",checks)
    write_csv(args.out/"development_pending_exclusions.csv",pending)
    write_csv(args.out/"development_sample_inventory.csv",sample_audit)

    with gzip.open(args.holdout,"rt",encoding="utf-8-sig",newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        rows = list(reader)
    assert header == ["Gene"]+[f"NL{i}" for i in range(1,7)]+[f"PSO{i}" for i in range(1,8)]
    assert len(set(r[0] for r in rows)) == len(rows) and all(r[0] for r in rows)
    values = np.array([[float(v) for v in r[1:]] for r in rows])
    assert np.isfinite(values).all() and (values>=0).all()
    libs = values.sum(axis=0)
    expr = np.log2(values/libs*1e6+1)
    geneidx = {r[0]:i for i,r in enumerate(rows)}
    holdchecks = []
    for oldr in holdold:
        ori = oldr["FrozenOrientation"]
        assert ori in ("UP_IN_PSORIASIS","DOWN_IN_PSORIASIS")
        sign = 1 if ori == "UP_IN_PSORIASIS" else -1
        x = sign*expr[geneidx[oldr["GeneSymbol"]]]
        m = metrics(x[6:],x[:6])
        delta = max(abs(m["AUC"]-float(oldr["DirectionalAUC"])),abs(m["AUC_SE"]-float(oldr["AUC_SE_DeLong"])),
                    abs(m["HedgesG"]-float(oldr["HedgesG_Aligned"])),abs(m["HedgesG_SE"]-float(oldr["HedgesG_SE"])))
        holdchecks.append(dict(PMID=oldr["PMID"], Gene=oldr["GeneSymbol"], **m,
                               PublishedDirectionConcordant=m["HedgesG"]>0,
                               PreHoldoutSignConcordant=m["HedgesG"]*float(oldr["PreHoldoutPooledHedgesG"])>0,
                               MaxAbsoluteDifference=delta, MatchesWithin1e_6=delta<1e-6))
    write_csv(args.out/"holdout_recomputed.csv",holdchecks)
    # Reproduce the observed alias collision as a diagnostic, never as the corrected estimate.
    if "MIF" in geneidx and "S100A9" in geneidx:
        combined = np.log2((values[geneidx["MIF"]]+values[geneidx["S100A9"]])/libs*1e6+1)
        collided = metrics(combined[6:],combined[:6])
        direct = metrics(expr[geneidx["S100A9"],6:],expr[geneidx["S100A9"],:6])
        source = next(r for r in holdold if r["GeneSymbol"] == "S100A9")
        write_csv(args.out/"S100A9_alias_collision_diagnostic.csv",[
            dict(Scenario="CURRENT_SYMBOL_S100A9_ONLY", **direct,
                 DifferenceFromUploadedG=direct["HedgesG"]-float(source["HedgesG_Aligned"])),
            dict(Scenario="DIAGNOSTIC_ONLY_S100A9_PLUS_DISTINCT_MIF_GENE", **collided,
                 DifferenceFromUploadedG=collided["HedgesG"]-float(source["HedgesG_Aligned"]))])
    write_csv(args.out/"holdout_column_diagnostics.csv",[dict(Column=header[j+1],
        AssignedGroup="HEALTHY_AS_IN_SUPPLIED_PROTOCOL" if j<6 else "PSORIASIS_AS_IN_SUPPLIED_PROTOCOL",
        Sum=libs[j], NonintegerValues=int((values[:,j]!=np.floor(values[:,j])).sum()),
        ZeroValues=int((values[:,j]==0).sum())) for j in range(13)])

    groups = defaultdict(list)
    for r in old:
        groups[(r["PMID"],r["GeneSymbol"])].append(r)
    fixed, meta_pending = [], []
    for p in pooled:
        key = (p["PMID"],p["GeneSymbol"])
        rs = groups[key]
        if {"GSE13355","GSE54456"} <= {r["GSE"] for r in rs} and excludes is None:
            meta_pending.append(dict(PMID=key[0], Gene=key[1], Reason="Additional synthesis-only overlap correction needs sample list"))
            continue
        synth = [dict(recomputed.get(r["TestID"], r)) for r in rs]
        both = {"GSE13355","GSE54456"} <= {r["GSE"] for r in rs}
        if both:
            for r in synth:
                if r["GSE"] != "GSE54456":
                    continue
                data = [x for x in matrices["GSE54456"] if x["GSM"] not in excludes]
                assert len(matrices["GSE54456"])-len(data) == 42
                sign = 1 if r["FrozenOrientation"] == "UP_IN_PSORIASIS" else -1
                mm = metrics([sign*float(x[r["GeneSymbol"]]) for x in data if x["FinalPhenotype"].startswith("PSORIASIS_LESIONAL_SKIN")],
                             [sign*float(x[r["GeneSymbol"]]) for x in data if x["FinalPhenotype"] == "HEALTHY_CONTROL_SKIN"])
                r.update(HedgesG_Aligned=mm["HedgesG"], HedgesG_SE=mm["HedgesG_SE"])
        fallback = sum(r["TestID"] not in recomputed for r in rs)
        m = corrected_meta(synth)
        fixed.append(dict(PMID=key[0], Gene=key[1], InputScope="PROCESSED_MATRICES_WITH_UPLOADED_ESTIMATE_FALLBACK",
                          UploadedEstimateFallbackN=fallback,
                          **m, OldTau2=p["Tau2_REML"], OldHedgesG=p["PooledHedgesG_REML"],
                          Tau2Change=m["Tau2"]-float(p["Tau2_REML"]),
                          HedgesGChange=m["HedgesG"]-float(p["PooledHedgesG_REML"]),
                          PrimaryClassification=classification(m["CI_Lower"],m["CI_Upper"])))
    write_csv(args.out/"reml_corrected_available_groups.csv",fixed)
    write_csv(args.out/"reml_pending_groups.csv",meta_pending)
    pis = []
    for p in pooled:
        k = int(p["K_IndependentCohorts"])
        mu = float(p["PooledHedgesG_REML"])
        tau = float(p["Tau2_REML"])
        crit = float(t.ppf(.975,k-1))
        seh = (float(p["mHK_CI95_Upper"])-float(p["mHK_CI95_Lower"]))/2/crit
        lo = mu-crit*math.sqrt(tau+seh*seh)
        hi = mu+crit*math.sqrt(tau+seh*seh)
        pis.append(dict(PMID=p["PMID"],Gene=p["GeneSymbol"],K=k,
                        Source="UPLOADED_STAGE4F_PARAMETERS_NOT_FULLY_REFITTED",
                        OldNormalPI_Lower=p["PredictionInterval_Lower"],OldNormalPI_Upper=p["PredictionInterval_Upper"],
                        NewPostHoc_mHK_t_PI_Lower=lo,NewPostHoc_mHK_t_PI_Upper=hi,
                        NewClassification=classification(lo,hi)))
    write_csv(args.out/"posthoc_prediction_interval_sensitivity.csv",pis)
    unique = {r["Gene"]:r for r in holdchecks}
    summary = dict(ReviewTimeUTC=datetime.now(timezone.utc).isoformat(), Python=platform.python_version(),
                   NumPy=np.__version__, SciPy=scipy.__version__, OriginalDevelopmentTests=len(old),
                   RecomputedDevelopmentTests=len(checks), PendingDevelopmentTests=len(pending),
                   DevelopmentMismatches=sum(not r["MatchesWithin1e_6"] for r in checks),
                   MaxDevelopmentDifference=max(r["MaxAbsoluteDifference"] for r in checks),
                   HoldoutClaims=len(holdchecks),HoldoutUniqueGenes=len(unique),
                   HoldoutMismatches=sum(not r["MatchesWithin1e_6"] for r in holdchecks),
                   MaxHoldoutDifference=max(r["MaxAbsoluteDifference"] for r in holdchecks),
                   HoldoutMedianAUC=statistics.median(r["AUC"] for r in holdchecks),
                   HoldoutMedianHedgesG=statistics.median(r["HedgesG"] for r in holdchecks),
                   HoldoutPublishedDirectionConcordant=sum(r["PublishedDirectionConcordant"] for r in unique.values()),
                   HoldoutPreSignConcordant=sum(r["PreHoldoutSignConcordant"] for r in unique.values()),
                   HoldoutUniquePerfectAUC=sum(r["AUC"]==1 for r in unique.values()),
                   HoldoutEffectCICrossesZero=[r["Gene"] for r in unique.values() if r["HedgesG_Lower"]<=0<=r["HedgesG_Upper"]],
                   REMLRecomputedFromSuppliedEstimates=len(fixed),REMLPending=len(meta_pending),
                   PostHocPIClassifications=dict(Counter(r["NewClassification"] for r in pis)),
                   HoldoutSHA256=hashlib.sha256(args.holdout.read_bytes()).hexdigest(),
                   ImportantScope="Historical original-grid arithmetic module. External reconstruction of target preprocessing and source eligibility is recorded separately under evidence/ and source_adjudicated_subset tables. This module alone does not validate source independence.")
    (args.out/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    paths = matrixfiles+[args.stage4a,args.stage4f,args.stage5b,args.holdout]
    if args.exclusions:
        paths.append(args.exclusions)
    write_csv(args.out/"input_checksums.csv",[dict(Filename=p.name,Bytes=p.stat().st_size,
               SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths])
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
