#requires -Version 7.0
$ErrorActionPreference='Stop'
Set-StrictMode -Version Latest
$Root=$env:PSORIASIS_CORRECTION_ROOT
if (-not $Root) { throw 'Set PSORIASIS_CORRECTION_ROOT to a separate correction workspace.' }
$ScriptDir=Join-Path $Root '04_scripts\stage4'
$LogDir=Join-Path $Root 'logs'
New-Item -ItemType Directory -Force -Path $ScriptDir,$LogDir | Out-Null
$Log=Join-Path $LogDir ("stage4d_auc_attenuation_random_effects_"+(Get-Date -Format 'yyyyMMdd_HHmmss')+".log")
Start-Transcript -Path $Log
Write-Host '=========================================================='
Write-Host 'STAGE 4D: AUC ATTENUATION + RANDOM-EFFECTS SYNTHESIS'
Write-Host '=========================================================='
Write-Host 'Uses only the Stage4D0 synthesis-independent frozen input.'
Write-Host 'Published AUC comparator mapping remains frozen from Stage4C2.'
Write-Host 'GSE295540 remains sealed.'
Write-Host ''
$Py=Join-Path $ScriptDir 'stage4d_auc_attenuation_random_effects.py'
@'

# CORRECTION PATCH 2026-09-07. Preserve original frozen results.

import csv, hashlib, math, statistics, zipfile
from collections import Counter, defaultdict
from pathlib import Path

import os
ROOT = Path(os.environ["PSORIASIS_CORRECTION_ROOT"])
S4D0=ROOT/"05_results"/"primary"/"stage4d0_synthesis_independence_freeze"
OUT=ROOT/"05_results"/"primary"/"stage4d_auc_attenuation_and_random_effects"

INPUT=S4D0/"single_gene_synthesis_ready_metrics_FINAL_FROZEN_STAGE4D0.csv"
COMPLETE0=S4D0/"STAGE4D0_COMPLETE.flag"

for p in [INPUT, COMPLETE0]:
    if not p.exists() or p.stat().st_size==0:
        raise RuntimeError(f"Required input missing/empty: {p}")

OUT.mkdir(parents=True,exist_ok=True)
COMPLETE=OUT/"STAGE4D_COMPLETE.flag"
if COMPLETE.exists():
    raise RuntimeError(f"Stage 4D is already complete/frozen: {COMPLETE}")

def clean(x): return (x or "").strip().strip('"').strip("'")
def rcsv(p):
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
def wcsv(p,rows,fields=None):
    rows=list(rows)
    if fields is None:
        fields=list(rows[0].keys()) if rows else []
    with p.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        if rows: w.writerows(rows)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):
            h.update(b)
    return h.hexdigest()

PROTOCOL="""PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT
STAGE 4D — PUBLISHED-AUC ATTENUATION + PAPER×GENE RANDOM-EFFECTS SYNTHESIS

Frozen input
Exactly:
05_results/primary/stage4d0_synthesis_independence_freeze/
single_gene_synthesis_ready_metrics_FINAL_FROZEN_STAGE4D0.csv

The Stage4D0 input has already enforced sample independence for synthesis where GSE13355
and GSE54456 would otherwise share biopsies.

PART A — PUBLISHED VS INDEPENDENT AUC
Published AUC comparators were frozen in Stage4C2 before attenuation was calculated.

Per validation row:
- AUC_Delta_IndependentMinusPublished = independent directional AUC - published comparator AUC.
- AUC_Attenuation_PublishedMinusIndependent = published comparator AUC - independent directional AUC.
  Positive values mean attenuation; negative values mean the independent AUC is higher.
- ExcessDiscriminationRetention = (independent AUC - 0.5)/(published AUC - 0.5).
  Values <0 are retained for directional reversals; values >1 mean stronger discrimination
  than the frozen published comparator.
- If published AUC = 0.5 exactly, ExcessDiscriminationRetention is undefined.

Attenuation is descriptive. No inferential variance is assigned to the published AUC because
the source claims do not all provide a common, transportable uncertainty estimate.

PART B — FORMAL CROSS-COHORT SYNTHESIS
Unit of synthesis: PMID × GeneSymbol.

Effect:
Aligned Hedges g from Stage4D0. Positive supports the prospectively frozen direction;
negative indicates reversal.

Variance:
Use the Stage4A/Stage4D0 HedgesG_SE squared.

Model:
Random-effects inverse-variance meta-analysis with between-cohort variance tau^2 estimated
by restricted maximum likelihood (REML).

REML objective, up to an additive constant:
0.5 * [sum(log(v_i + tau^2)) + log(sum(w_i)) + sum(w_i*(y_i-mu)^2)]
where w_i = 1/(v_i + tau^2) and mu = sum(w_i*y_i)/sum(w_i).

Tau^2 is constrained to >=0 and optimized numerically.
Pooled 95% CI uses normal critical value 1.96:
mu +/- 1.96*sqrt(1/sum(w_i)).
Prediction interval for k>=3:
mu +/- 1.96*sqrt(tau^2 + 1/sum(w_i)).

Heterogeneity:
- Cochran Q from fixed-effect inverse-variance weights.
- I^2 = max(0,(Q-df)/Q)*100 when Q>0.
- tau^2 from REML.

Synthesis interpretation:
- REPLICATED_CI_EXCLUDES_ZERO: pooled g >0 and lower 95% CI >0.
- POSITIVE_INCONCLUSIVE: pooled g >0 but CI includes 0.
- REVERSED_CI_EXCLUDES_ZERO: pooled g <0 and upper 95% CI <0.
- NEGATIVE_INCONCLUSIVE: pooled g <0 but CI includes 0.
- NULL_OR_TIE: pooled g =0.

No marker is removed for heterogeneity or reversal.
No global meta-analysis across all 39 markers is performed because marker-level syntheses
reuse many of the same cohorts and are not mutually independent.

No threshold optimization, post-hoc sign flipping, feature selection, or GSE295540 access.
"""
(OUT/"STAGE4D_SYNTHESIS_PROTOCOL_FROZEN.txt").write_text(PROTOCOL,encoding="utf-8")

rows=rcsv(INPUT)
if len(rows)!=202:
    raise RuntimeError(f"Expected 202 synthesis-ready rows; found {len(rows)}")
if any(clean(r.get("GSE295540Accessed",""))=="YES" for r in rows):
    raise RuntimeError("GSE295540 access flag detected.")
if any(clean(r.get("SynthesisSampleIndependentFromOtherIncludedBenchmarks",""))!="YES" for r in rows):
    raise RuntimeError("At least one Stage4D0 row is not synthesis-independent.")
if any(clean(r.get("PublishedComparatorAUC",""))=="" for r in rows):
    raise RuntimeError("At least one row lacks frozen published comparator.")

# ---------- REML implementation ----------
def reml_components(y,v,tau2):
    w=[1.0/(vi+tau2) for vi in v]
    sw=sum(w)
    mu=sum(wi*yi for wi,yi in zip(w,y))/sw
    q=sum(wi*(yi-mu)**2 for wi,yi in zip(w,y))
    obj=0.5*(sum(math.log(vi+tau2) for vi in v)+math.log(sw)+q)
    return obj,mu,w,q

def golden_minimize(f,a,b,tol=1e-12,maxiter=300):
    gr=(math.sqrt(5)-1)/2
    c=b-gr*(b-a)
    d=a+gr*(b-a)
    fc=f(c); fd=f(d)
    for _ in range(maxiter):
        if abs(b-a) <= tol*(1+abs(a)+abs(b)):
            break
        if fc < fd:
            b=d; d=c; fd=fc
            c=b-gr*(b-a); fc=f(c)
        else:
            a=c; c=d; fc=fd
            d=a+gr*(b-a); fd=f(d)
    x=(a+b)/2
    return x,f(x)

def estimate_tau2_reml(y,v):
    k=len(y)
    if k<2:
        return 0.0
    f=lambda t: reml_components(y,v,t)[0]
    f0=f(0.0)

    # Build an adaptive upper bound well beyond the observed effect spread.
    var_y=statistics.variance(y) if k>=2 else 0.0
    mean_v=statistics.mean(v)
    upper=max(1e-8, var_y, mean_v, 0.01)
    prev=f(upper)
    for _ in range(30):
        nxt=upper*4.0
        fn=f(nxt)
        if fn >= prev:
            upper = nxt  # correction: optimum can lie between upper and nxt
            break
        upper=nxt
        prev=fn

    t_hat,obj=golden_minimize(f,0.0,upper)
    # Boundary REML solution at zero is allowed.
    if f0 <= obj + 1e-10:
        return 0.0
    return max(0.0,t_hat)

def meta_reml(y,se):
    k=len(y)
    if k<1:
        raise RuntimeError("Empty meta-analysis group")
    v=[s*s for s in se]
    if any((not math.isfinite(x) or x<=0) for x in v):
        raise RuntimeError(f"Invalid sampling variance(s): {v}")

    if k==1:
        mu=y[0]; sem=se[0]
        return {
            "k":1,"tau2":0.0,"Q":0.0,"df":0,"I2":0.0,
            "mu":mu,"se":sem,"lo":mu-1.96*sem,"hi":mu+1.96*sem,
            "pred_lo":"","pred_hi":""
        }

    # Fixed-effect Q for heterogeneity.
    wf=[1.0/vi for vi in v]
    swf=sum(wf)
    muf=sum(wi*yi for wi,yi in zip(wf,y))/swf
    Q=sum(wi*(yi-muf)**2 for wi,yi in zip(wf,y))
    df=k-1
    I2=max(0.0,(Q-df)/Q*100.0) if Q>0 else 0.0

    tau2=estimate_tau2_reml(y,v)
    _,mu,w,_=reml_components(y,v,tau2)
    sw=sum(w)
    sem=math.sqrt(1.0/sw)
    lo=mu-1.96*sem
    hi=mu+1.96*sem
    if k>=3:
        predse=math.sqrt(tau2+sem*sem)
        plo=mu-1.96*predse
        phi=mu+1.96*predse
    else:
        plo=phi=""
    return {
        "k":k,"tau2":tau2,"Q":Q,"df":df,"I2":I2,
        "mu":mu,"se":sem,"lo":lo,"hi":hi,
        "pred_lo":plo,"pred_hi":phi
    }

# ---------- Part A: row-level attenuation ----------
atten=[]
for r in rows:
    ind=float(r["AUC_Directional"])
    pub=float(r["PublishedComparatorAUC"])
    delta=ind-pub
    att=pub-ind
    if abs(pub-0.5)<1e-15:
        retention=""
    else:
        retention=(ind-0.5)/(pub-0.5)
    if att>0.05:
        band="ATTENUATION_GT_0.05"
    elif att<-0.05:
        band="INDEPENDENT_HIGHER_GT_0.05"
    else:
        band="WITHIN_0.05"
    atten.append({
        "TestID":r["TestID"],
        "PMID":r["PMID"],
        "GeneSymbol":r["GeneSymbol"],
        "GSE":r["GSE"],
        "BenchmarkPlatformClass":clean(r.get("BenchmarkPlatformClass","")),
        "FrozenOrientation":r["FrozenOrientation"],
        "PublishedComparatorAUC":f"{pub:.12g}",
        "IndependentDirectionalAUC":f"{ind:.12g}",
        "AUC_Delta_IndependentMinusPublished":f"{delta:.12g}",
        "AUC_Attenuation_PublishedMinusIndependent":f"{att:.12g}",
        "ExcessDiscriminationRetention":"" if retention=="" else f"{retention:.12g}",
        "AttenuationBand":band,
        "SynthesisMetricSource":r["SynthesisMetricSource"],
    })
wcsv(OUT/"published_vs_independent_auc_attenuation_STAGE4D.csv",atten)

# ---------- Part B: paper×gene summaries and meta-analysis ----------
groups=defaultdict(list)
for r in rows:
    groups[(clean(r["PMID"]),clean(r["GeneSymbol"]))].append(r)
if len(groups)!=39:
    raise RuntimeError(f"Expected 39 paper×gene groups; found {len(groups)}")

meta=[]
pg_atten=[]
for (pmid,gene),rs in sorted(groups.items()):
    ys=[float(r["HedgesG_Aligned"]) for r in rs]
    ses=[float(r["HedgesG_SE"]) for r in rs]
    m=meta_reml(ys,ses)

    if m["mu"]>0 and m["lo"]>0:
        synth="REPLICATED_CI_EXCLUDES_ZERO"
    elif m["mu"]>0:
        synth="POSITIVE_INCONCLUSIVE"
    elif m["mu"]<0 and m["hi"]<0:
        synth="REVERSED_CI_EXCLUDES_ZERO"
    elif m["mu"]<0:
        synth="NEGATIVE_INCONCLUSIVE"
    else:
        synth="NULL_OR_TIE"

    aucs=[float(r["AUC_Directional"]) for r in rs]
    pubs=[float(r["PublishedComparatorAUC"]) for r in rs]
    atts=[p-a for p,a in zip(pubs,aucs)]
    deltas=[a-p for p,a in zip(pubs,aucs)]
    rets=[(a-0.5)/(p-0.5) for p,a in zip(pubs,aucs) if abs(p-0.5)>1e-15]

    meta.append({
        "PMID":pmid,
        "GeneSymbol":gene,
        "K_IndependentCohorts":m["k"],
        "Benchmarks":"|".join(sorted(r["GSE"] for r in rs)),
        "PooledHedgesG_REML":f"{m['mu']:.12g}",
        "PooledHedgesG_SE":f"{m['se']:.12g}",
        "PooledHedgesG_CI95_Lower":f"{m['lo']:.12g}",
        "PooledHedgesG_CI95_Upper":f"{m['hi']:.12g}",
        "Tau2_REML":f"{m['tau2']:.12g}",
        "CochranQ":f"{m['Q']:.12g}",
        "Q_df":m["df"],
        "I2_Percent":f"{m['I2']:.12g}",
        "PredictionInterval_Lower":"" if m["pred_lo"]=="" else f"{m['pred_lo']:.12g}",
        "PredictionInterval_Upper":"" if m["pred_hi"]=="" else f"{m['pred_hi']:.12g}",
        "SynthesisClassification":synth,
        "MedianIndependentDirectionalAUC":f"{statistics.median(aucs):.12g}",
        "MinIndependentDirectionalAUC":f"{min(aucs):.12g}",
        "MaxIndependentDirectionalAUC":f"{max(aucs):.12g}",
        "DirectionalAUC_Below_0.5_Count":sum(a<0.5 for a in aucs),
    })

    pg_atten.append({
        "PMID":pmid,
        "GeneSymbol":gene,
        "K_IndependentCohorts":len(rs),
        "PublishedComparatorAUCs":"|".join(f"{x:.12g}" for x in sorted(set(pubs))),
        "MedianPublishedComparatorAUC":f"{statistics.median(pubs):.12g}",
        "MedianIndependentDirectionalAUC":f"{statistics.median(aucs):.12g}",
        "MedianAUC_Delta_IndependentMinusPublished":f"{statistics.median(deltas):.12g}",
        "MedianAUC_Attenuation_PublishedMinusIndependent":f"{statistics.median(atts):.12g}",
        "MinAUC_Attenuation":f"{min(atts):.12g}",
        "MaxAUC_Attenuation":f"{max(atts):.12g}",
        "MedianExcessDiscriminationRetention":"" if not rets else f"{statistics.median(rets):.12g}",
        "CohortsWithAttenuation_GT_0.05":sum(x>0.05 for x in atts),
        "CohortsWithin_0.05":sum(abs(x)<=0.05 for x in atts),
        "CohortsIndependentHigher_GT_0.05":sum(x<-0.05 for x in atts),
    })

wcsv(OUT/"paper_gene_random_effects_hedges_g_STAGE4D.csv",meta)
wcsv(OUT/"paper_gene_auc_attenuation_summary_STAGE4D.csv",pg_atten)

# ---------- Descriptive overall summaries ----------
row_att=[float(r["AUC_Attenuation_PublishedMinusIndependent"]) for r in atten]
row_delta=[float(r["AUC_Delta_IndependentMinusPublished"]) for r in atten]
pg_med_att=[float(r["MedianAUC_Attenuation_PublishedMinusIndependent"]) for r in pg_atten]
pooled_g=[float(r["PooledHedgesG_REML"]) for r in meta]
i2=[float(r["I2_Percent"]) for r in meta]
class_counts=Counter(r["SynthesisClassification"] for r in meta)
band_counts=Counter(r["AttenuationBand"] for r in atten)

# Rank useful tables.
largest_att=sorted(atten,key=lambda r:float(r["AUC_Attenuation_PublishedMinusIndependent"]),reverse=True)[:25]
wcsv(OUT/"top25_largest_auc_attenuations_STAGE4D.csv",largest_att)

lowest_meta=sorted(meta,key=lambda r:float(r["PooledHedgesG_REML"]))
wcsv(OUT/"paper_gene_synthesis_ranked_by_pooled_g_STAGE4D.csv",lowest_meta)

summary=[
"PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT",
"STAGE 4D COMPLETE — AUC ATTENUATION + RANDOM-EFFECTS SYNTHESIS","",
"Input Stage4D0 synthesis-independent rows: 202 / 202",
"Paper×gene syntheses: 39 / 39",
"GSE295540 accessed: NO",
"Post-hoc sign flipping: NO",
"Performance-driven exclusions: NO",
"Grand meta-analysis across markers: NO","",
"AUC ATTENUATION — ROW LEVEL",
f"Median published-minus-independent AUC attenuation: {statistics.median(row_att):.6f}",
f"Mean published-minus-independent AUC attenuation: {statistics.mean(row_att):.6f}",
f"Median independent-minus-published AUC delta: {statistics.median(row_delta):.6f}",
f"Rows with attenuation >0.05: {band_counts['ATTENUATION_GT_0.05']} / 202",
f"Rows within +/-0.05: {band_counts['WITHIN_0.05']} / 202",
f"Rows with independent AUC > published by >0.05: {band_counts['INDEPENDENT_HIGHER_GT_0.05']} / 202","",
"AUC ATTENUATION — PAPER×GENE WEIGHTED DESCRIPTIVE SUMMARY",
f"Median of the 39 paper×gene median attenuations: {statistics.median(pg_med_att):.6f}",
f"Paper×gene median attenuation range: {min(pg_med_att):.6f} to {max(pg_med_att):.6f}","",
"RANDOM-EFFECTS HEDGES g SYNTHESIS",
]
for k in ["REPLICATED_CI_EXCLUDES_ZERO","POSITIVE_INCONCLUSIVE",
          "REVERSED_CI_EXCLUDES_ZERO","NEGATIVE_INCONCLUSIVE","NULL_OR_TIE"]:
    summary.append(f"{k}: {class_counts[k]}")
summary += [
f"Median pooled aligned Hedges g across 39 marker syntheses: {statistics.median(pooled_g):.6f}",
f"Pooled g range: {min(pooled_g):.6f} to {max(pooled_g):.6f}",
f"Median I2 across marker syntheses: {statistics.median(i2):.2f}%",
"",
"Lowest pooled-g syntheses:"
]
for r in lowest_meta[:8]:
    summary.append(
        f"  PMID {r['PMID']} / {r['GeneSymbol']}: pooled g={float(r['PooledHedgesG_REML']):.6f} "
        f"[{float(r['PooledHedgesG_CI95_Lower']):.6f}, {float(r['PooledHedgesG_CI95_Upper']):.6f}], "
        f"I2={float(r['I2_Percent']):.1f}%, class={r['SynthesisClassification']}"
    )
summary += [
"",
"INTERPRETIVE LIMIT",
"Published-vs-independent AUC attenuation is descriptive because published AUC uncertainty",
"is not harmonized across claims. Formal inference is therefore based on cross-cohort aligned",
"Hedges g, not on treating the frozen published AUC as an error-free meta-analytic observation.",
"",
"NEXT GATE",
"Stage 4E can produce manuscript-ready primary tables/figures and prespecified sensitivity",
"analyses while leaving GSE295540 sealed for the final temporal holdout."
]
(OUT/"STAGE4D_SUMMARY.txt").write_text("\n".join(summary)+"\n",encoding="utf-8")

manifest=[]
for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {"STAGE4D_FILE_MANIFEST.csv","Stage4D_AUCAttenuation_RandomEffects_Bundle.zip","STAGE4D_COMPLETE.flag"}:
        manifest.append({"RelativePath":str(p.relative_to(OUT)),"Bytes":p.stat().st_size,"SHA256":sha(p)})
wcsv(OUT/"STAGE4D_FILE_MANIFEST.csv",manifest)

bundle=OUT/"Stage4D_AUCAttenuation_RandomEffects_Bundle.zip"
if bundle.exists(): bundle.unlink()
with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p!=bundle:
            z.write(p,p.relative_to(OUT))

COMPLETE.write_text("Completed. 202 attenuation rows + 39 REML paper-gene syntheses. GSE295540 not accessed.\n",encoding="utf-8")

print("\n==========================================================")
print("STAGE 4D COMPLETE")
print("==========================================================")
print("\n".join(summary))
print("\nResults:")
print(OUT/"published_vs_independent_auc_attenuation_STAGE4D.csv")
print(OUT/"paper_gene_random_effects_hedges_g_STAGE4D.csv")
print("Bundle:")
print(bundle)


'@ | Set-Content -Encoding UTF8 -Path $Py
$PyExe=$null
if(Get-Command py.exe -ErrorAction SilentlyContinue){$PyExe='py.exe'}
elseif(Get-Command python.exe -ErrorAction SilentlyContinue){$PyExe='python.exe'}
elseif(Get-Command python -ErrorAction SilentlyContinue){$PyExe='python'}
else{throw 'Python 3 was not found.'}
if($PyExe -eq 'py.exe'){
    & $PyExe -3 -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 4D syntax check failed.'}
    & $PyExe -3 $Py
}else{
    & $PyExe -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 4D syntax check failed.'}
    & $PyExe $Py
}
if($LASTEXITCODE -ne 0){throw "Stage 4D failed with exit code $LASTEXITCODE"}
Stop-Transcript
