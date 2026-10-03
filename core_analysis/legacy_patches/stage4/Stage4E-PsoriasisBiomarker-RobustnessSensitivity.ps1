#requires -Version 7.0
$ErrorActionPreference='Stop'
Set-StrictMode -Version Latest
$Root=$env:PSORIASIS_CORRECTION_ROOT
if (-not $Root) { throw 'Set PSORIASIS_CORRECTION_ROOT to a separate correction workspace.' }
$ScriptDir=Join-Path $Root '04_scripts\stage4'
$LogDir=Join-Path $Root 'logs'
New-Item -ItemType Directory -Force -Path $ScriptDir,$LogDir | Out-Null
$Log=Join-Path $LogDir ("stage4e_primary_robustness_"+(Get-Date -Format 'yyyyMMdd_HHmmss')+".log")
Start-Transcript -Path $Log
Write-Host '=========================================================='
Write-Host 'STAGE 4E: PRIMARY ROBUSTNESS / SENSITIVITY ANALYSES'
Write-Host '=========================================================='
Write-Host 'Primary Stage4D results remain unchanged.'
Write-Host 'LOCO + platform strata + small-control exclusion + mHK + prediction intervals.'
Write-Host 'GSE295540 remains sealed.'
Write-Host ''
$Py=Join-Path $ScriptDir 'stage4e_primary_robustness.py'
@'

# CORRECTION PATCH 2026-09-07. Preserve original frozen results.

import csv, hashlib, math, statistics, zipfile
from collections import Counter, defaultdict
from pathlib import Path

import os
ROOT = Path(os.environ["PSORIASIS_CORRECTION_ROOT"])
S4D0=ROOT/"05_results"/"primary"/"stage4d0_synthesis_independence_freeze"
S4D=ROOT/"05_results"/"primary"/"stage4d_auc_attenuation_and_random_effects"
OUT=ROOT/"05_results"/"sensitivity"/"stage4e_primary_robustness"

INPUT=S4D0/"single_gene_synthesis_ready_metrics_FINAL_FROZEN_STAGE4D0.csv"
PRIMARY=S4D/"paper_gene_random_effects_hedges_g_STAGE4D.csv"
ATTEN=S4D/"published_vs_independent_auc_attenuation_STAGE4D.csv"

for p in [INPUT,PRIMARY,ATTEN,S4D0/"STAGE4D0_COMPLETE.flag",S4D/"STAGE4D_COMPLETE.flag"]:
    if not p.exists() or p.stat().st_size==0:
        raise RuntimeError(f"Required input missing/empty: {p}")

OUT.mkdir(parents=True,exist_ok=True)
COMPLETE=OUT/"STAGE4E_COMPLETE.flag"
if COMPLETE.exists():
    raise RuntimeError(f"Stage 4E already complete/frozen: {COMPLETE}")

def clean(x): return (x or "").strip().strip('"').strip("'")
def rcsv(p):
    with p.open(newline="",encoding="utf-8-sig") as f: return list(csv.DictReader(f))
def wcsv(p,rows,fields=None):
    rows=list(rows)
    if fields is None: fields=list(rows[0].keys()) if rows else []
    with p.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        if rows: w.writerows(rows)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

PROTOCOL="""PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT
STAGE 4E — PRIMARY SYNTHESIS ROBUSTNESS / SENSITIVITY ANALYSES

Primary Stage4D results remain unchanged.

Input:
Stage4D0 synthesis-independent per-cohort metrics.

Sensitivity A — Leave-one-cohort-out (LOCO)
For every PMID × GeneSymbol with >=3 cohorts, repeat the same REML Hedges-g synthesis after
removing each cohort in turn. No cohort is removed because of its performance.
Report whether the pooled sign and normal-theory 95% CI classification remain stable.

Sensitivity B — Platform-stratified synthesis
Frozen platform classes:
MICROARRAY = GSE13355, GSE14905, GSE78097, GSE201827
RNASEQ = GSE54456, GSE66511, GSE121212
For strata with >=2 cohorts, estimate separate REML pooled Hedges g.
This is descriptive/robustness analysis; no between-platform hypothesis test is performed.

Sensitivity C — Predefined small-control exclusion
Stage4B predefined SMALL_HEALTHY_GROUP as N_Healthy < 10.
Repeat each paper×gene synthesis after excluding such rows.
In this benchmark set this is expected to remove GSE78097 (6 healthy controls).
This sensitivity is based on sample size, not performance.

Sensitivity D — Modified Hartung-Knapp (mHK) interval
Keep the primary REML tau^2 and random-effects weights.
For k>=2:
q = sum(w_i*(y_i-mu)^2)/(k-1)
q* = max(1, q)
SE_mHK = sqrt(q*/sum(w_i))
CI = mu +/- t_(k-1,0.975)*SE_mHK.
The max(1,q) modification prevents the mHK interval from becoming narrower than the
ordinary random-effects interval due to q<1.
For the available maximum of seven cohorts, exact two-sided 95% t critical values for df 1..6
are frozen in the script.

Sensitivity E — Prediction-interval interpretation
Use the Stage4D REML prediction intervals. Classify each paper×gene as:
- PI_POSITIVE if lower PI > 0
- PI_CROSSES_ZERO if lower <=0<= upper
- PI_NEGATIVE if upper PI < 0
Prediction interval is unavailable for k<3.

AUC attenuation sensitivity
Report the Stage4D row-level attenuation distribution after excluding SMALL_HEALTHY_GROUP rows.
No new published comparator selection is performed.

No primary results are removed or replaced.
No post-hoc sign flipping.
No GSE295540 access.
"""
(OUT/"STAGE4E_SENSITIVITY_PROTOCOL_FROZEN.txt").write_text(PROTOCOL,encoding="utf-8")

rows=rcsv(INPUT)
primary=rcsv(PRIMARY)
atten=rcsv(ATTEN)
if len(rows)!=202 or len(primary)!=39 or len(atten)!=202:
    raise RuntimeError(f"Unexpected input sizes: rows={len(rows)} primary={len(primary)} attenuation={len(atten)}")

MICRO={"GSE13355","GSE14905","GSE78097","GSE201827"}
RNA={"GSE54456","GSE66511","GSE121212"}

def platform(g):
    if g in MICRO: return "MICROARRAY"
    if g in RNA: return "RNASEQ"
    raise RuntimeError(f"Unknown platform for {g}")

def reml_components(y,v,tau2):
    w=[1/(vi+tau2) for vi in v]
    sw=sum(w)
    mu=sum(wi*yi for wi,yi in zip(w,y))/sw
    q=sum(wi*(yi-mu)**2 for wi,yi in zip(w,y))
    obj=.5*(sum(math.log(vi+tau2) for vi in v)+math.log(sw)+q)
    return obj,mu,w,q

def golden(f,a,b,tol=1e-12,maxiter=300):
    gr=(math.sqrt(5)-1)/2
    c=b-gr*(b-a); d=a+gr*(b-a); fc=f(c); fd=f(d)
    for _ in range(maxiter):
        if abs(b-a)<=tol*(1+abs(a)+abs(b)): break
        if fc<fd:
            b=d; d=c; fd=fc; c=b-gr*(b-a); fc=f(c)
        else:
            a=c; c=d; fc=fd; d=a+gr*(b-a); fd=f(d)
    x=(a+b)/2
    return x,f(x)

def tau_reml(y,v):
    if len(y)<2: return 0.0
    f=lambda t:reml_components(y,v,t)[0]
    f0=f(0.0)
    vy=statistics.variance(y)
    upper=max(1e-8,vy,statistics.mean(v),.01)
    prev=f(upper)
    for _ in range(30):
        nxt=upper*4
        fn=f(nxt)
        if fn>=prev:
            upper = nxt  # correction: retain the rising endpoint in the bracket
            break
        upper=nxt; prev=fn
    th,obj=golden(f,0,upper)
    return 0.0 if f0<=obj+1e-10 else max(0.0,th)

def meta(rs):
    y=[float(r["HedgesG_Aligned"]) for r in rs]
    se=[float(r["HedgesG_SE"]) for r in rs]
    v=[s*s for s in se]
    k=len(y)
    if k==1:
        mu=y[0]; sem=se[0]
        return {"k":1,"mu":mu,"se":sem,"lo":mu-1.96*sem,"hi":mu+1.96*sem,
                "tau2":0.0,"I2":0.0,"Q":0.0,"pred_lo":"","pred_hi":"","weights":[1/v[0]],"q_re":0.0}
    wf=[1/x for x in v]; sw=sum(wf); muf=sum(wi*yi for wi,yi in zip(wf,y))/sw
    Q=sum(wi*(yi-muf)**2 for wi,yi in zip(wf,y)); df=k-1
    I2=max(0,(Q-df)/Q*100) if Q>0 else 0
    t2=tau_reml(y,v)
    _,mu,w,qre=reml_components(y,v,t2)
    sem=math.sqrt(1/sum(w))
    lo=mu-1.96*sem; hi=mu+1.96*sem
    if k>=3:
        pse=math.sqrt(t2+sem*sem); plo=mu-1.96*pse; phi=mu+1.96*pse
    else: plo=phi=""
    return {"k":k,"mu":mu,"se":sem,"lo":lo,"hi":hi,"tau2":t2,"I2":I2,"Q":Q,
            "pred_lo":plo,"pred_hi":phi,"weights":w,"q_re":qre}

def klass(m):
    if m["mu"]>0 and m["lo"]>0: return "REPLICATED_CI_EXCLUDES_ZERO"
    if m["mu"]>0: return "POSITIVE_INCONCLUSIVE"
    if m["mu"]<0 and m["hi"]<0: return "REVERSED_CI_EXCLUDES_ZERO"
    if m["mu"]<0: return "NEGATIVE_INCONCLUSIVE"
    return "NULL_OR_TIE"

T975={1:12.7062047364,2:4.30265272991,3:3.18244630528,4:2.7764451052,5:2.57058183564,6:2.44691185114}

groups=defaultdict(list)
for r in rows:
    groups[(clean(r["PMID"]),clean(r["GeneSymbol"]))].append(r)
if len(groups)!=39: raise RuntimeError(f"Expected 39 groups, found {len(groups)}")

pmap={(r["PMID"],r["GeneSymbol"]):r for r in primary}

# A: LOCO
loco=[]
loco_summary=[]
for key,rs in sorted(groups.items()):
    base=meta(rs)
    runs=[]
    if len(rs)>=3:
        for drop in sorted(r["GSE"] for r in rs):
            kept=[r for r in rs if r["GSE"]!=drop]
            m=meta(kept)
            rr={
                "PMID":key[0],"GeneSymbol":key[1],"DroppedGSE":drop,"K_Remaining":m["k"],
                "PooledHedgesG_REML":f"{m['mu']:.12g}",
                "CI95_Lower":f"{m['lo']:.12g}","CI95_Upper":f"{m['hi']:.12g}",
                "Tau2_REML":f"{m['tau2']:.12g}","I2_Percent":f"{m['I2']:.12g}",
                "Classification":klass(m)
            }
            loco.append(rr); runs.append(rr)
    if runs:
        classes=Counter(r["Classification"] for r in runs)
        all_pos=all(float(r["PooledHedgesG_REML"])>0 for r in runs)
        all_pos_sig=all(r["Classification"]=="REPLICATED_CI_EXCLUDES_ZERO" for r in runs)
        all_neg_sig=all(r["Classification"]=="REVERSED_CI_EXCLUDES_ZERO" for r in runs)
        min_g=min(float(r["PooledHedgesG_REML"]) for r in runs)
        max_g=max(float(r["PooledHedgesG_REML"]) for r in runs)
        loco_summary.append({
            "PMID":key[0],"GeneSymbol":key[1],"K_Full":len(rs),"LOCO_Runs":len(runs),
            "FullClassification":klass(base),
            "AllLOCO_PooledSignsPositive":"YES" if all_pos else "NO",
            "AllLOCO_Replicated_CI_Excludes_Zero":"YES" if all_pos_sig else "NO",
            "AllLOCO_Reversed_CI_Excludes_Zero":"YES" if all_neg_sig else "NO",
            "MinLOCO_PooledG":f"{min_g:.12g}","MaxLOCO_PooledG":f"{max_g:.12g}",
            "LOCO_Classifications":"|".join(f"{k}:{v}" for k,v in sorted(classes.items()))
        })
wcsv(OUT/"leave_one_cohort_out_results_STAGE4E.csv",loco)
wcsv(OUT/"leave_one_cohort_out_summary_STAGE4E.csv",loco_summary)

# B: platform strata
plat_rows=[]
plat_summary=[]
for key,rs in sorted(groups.items()):
    rec={"PMID":key[0],"GeneSymbol":key[1]}
    signs={}
    for pname,pool in [("MICROARRAY",MICRO),("RNASEQ",RNA)]:
        sub=[r for r in rs if r["GSE"] in pool]
        if len(sub)>=2:
            m=meta(sub)
            plat_rows.append({
                "PMID":key[0],"GeneSymbol":key[1],"Platform":pname,"K":len(sub),
                "Benchmarks":"|".join(sorted(r["GSE"] for r in sub)),
                "PooledHedgesG_REML":f"{m['mu']:.12g}",
                "CI95_Lower":f"{m['lo']:.12g}","CI95_Upper":f"{m['hi']:.12g}",
                "Tau2_REML":f"{m['tau2']:.12g}","I2_Percent":f"{m['I2']:.12g}",
                "Classification":klass(m)
            })
            signs[pname]=1 if m["mu"]>0 else -1 if m["mu"]<0 else 0
    if signs:
        plat_summary.append({
            "PMID":key[0],"GeneSymbol":key[1],
            "MicroarrayAvailable":"YES" if "MICROARRAY" in signs else "NO",
            "RNAseqAvailable":"YES" if "RNASEQ" in signs else "NO",
            "PlatformPooledSignsConcordant":"YES" if len(signs)==2 and len(set(signs.values()))==1 else
                                           "NO" if len(signs)==2 else "NOT_BOTH_AVAILABLE"
        })
wcsv(OUT/"platform_stratified_random_effects_STAGE4E.csv",plat_rows)
wcsv(OUT/"platform_sign_concordance_STAGE4E.csv",plat_summary)

# C: exclude N_Healthy <10
small=[]
small_summary=[]
for key,rs in sorted(groups.items()):
    kept=[r for r in rs if int(r["N_Healthy"])>=10]
    dropped=[r for r in rs if int(r["N_Healthy"])<10]
    if not dropped:
        continue
    if not kept:
        raise RuntimeError(f"{key}: small-control sensitivity removes all rows")
    m=meta(kept)
    small.append({
        "PMID":key[0],"GeneSymbol":key[1],"K_Original":len(rs),"K_After":len(kept),
        "DroppedBenchmarks":"|".join(sorted(r["GSE"] for r in dropped)),
        "PooledHedgesG_REML":f"{m['mu']:.12g}",
        "CI95_Lower":f"{m['lo']:.12g}","CI95_Upper":f"{m['hi']:.12g}",
        "Tau2_REML":f"{m['tau2']:.12g}","I2_Percent":f"{m['I2']:.12g}",
        "Classification":klass(m)
    })
wcsv(OUT/"small_healthy_group_exclusion_sensitivity_STAGE4E.csv",small)

# D: modified Hartung-Knapp
mhk=[]
for key,rs in sorted(groups.items()):
    m=meta(rs); k=m["k"]
    if k<2: continue
    df=k-1
    if df not in T975: raise RuntimeError(f"No frozen t critical for df={df}")
    qstar=max(1.0,m["q_re"]/df)
    se_hk=math.sqrt(qstar/sum(m["weights"]))
    crit=T975[df]
    lo=m["mu"]-crit*se_hk; hi=m["mu"]+crit*se_hk
    mm={"mu":m["mu"],"lo":lo,"hi":hi}
    c=klass(mm)
    mhk.append({
        "PMID":key[0],"GeneSymbol":key[1],"K":k,
        "PooledHedgesG_REML":f"{m['mu']:.12g}",
        "Tau2_REML":f"{m['tau2']:.12g}",
        "mHK_qstar":f"{qstar:.12g}","mHK_SE":f"{se_hk:.12g}",
        "t_Critical_975":f"{crit:.12g}",
        "mHK_CI95_Lower":f"{lo:.12g}","mHK_CI95_Upper":f"{hi:.12g}",
        "PrimaryNormalCIClassification":klass(m),
        "mHK_Classification":c
    })
wcsv(OUT/"modified_hartung_knapp_sensitivity_STAGE4E.csv",mhk)

# E: prediction intervals from primary Stage4D
pi=[]
for r in primary:
    plo=clean(r["PredictionInterval_Lower"]); phi=clean(r["PredictionInterval_Upper"])
    if plo=="" or phi=="":
        pc="UNAVAILABLE_K_LT_3"
    else:
        a=float(plo); b=float(phi)
        if a>0: pc="PI_POSITIVE"
        elif b<0: pc="PI_NEGATIVE"
        else: pc="PI_CROSSES_ZERO"
    pi.append({
        "PMID":r["PMID"],"GeneSymbol":r["GeneSymbol"],"K_IndependentCohorts":r["K_IndependentCohorts"],
        "PooledHedgesG_REML":r["PooledHedgesG_REML"],
        "PredictionInterval_Lower":plo,"PredictionInterval_Upper":phi,
        "PredictionIntervalClassification":pc
    })
wcsv(OUT/"prediction_interval_interpretation_STAGE4E.csv",pi)

# AUC attenuation excluding small healthy rows.
small_ids={clean(r["TestID"]) for r in rows if int(r["N_Healthy"])<10}
att_kept=[r for r in atten if clean(r["TestID"]) not in small_ids]
atts=[float(r["AUC_Attenuation_PublishedMinusIndependent"]) for r in att_kept]
bands=Counter(
    "ATTENUATION_GT_0.05" if x>0.05 else "INDEPENDENT_HIGHER_GT_0.05" if x<-0.05 else "WITHIN_0.05"
    for x in atts
)
wcsv(OUT/"auc_attenuation_small_control_exclusion_summary_STAGE4E.csv",[{
    "RowsOriginal":len(atten),"RowsAfterExcluding_NHealthy_lt10":len(att_kept),
    "ExcludedRows":len(atten)-len(att_kept),
    "MedianPublishedMinusIndependentAttenuation":f"{statistics.median(atts):.12g}",
    "MeanPublishedMinusIndependentAttenuation":f"{statistics.mean(atts):.12g}",
    "Attenuation_GT_0.05":bands["ATTENUATION_GT_0.05"],
    "Within_0.05":bands["WITHIN_0.05"],
    "IndependentHigher_GT_0.05":bands["INDEPENDENT_HIGHER_GT_0.05"]
}])

# Summary
loco_all_rep=sum(r["AllLOCO_Replicated_CI_Excludes_Zero"]=="YES" for r in loco_summary)
loco_all_rev=sum(r["AllLOCO_Reversed_CI_Excludes_Zero"]=="YES" for r in loco_summary)
mhk_counts=Counter(r["mHK_Classification"] for r in mhk)
pi_counts=Counter(r["PredictionIntervalClassification"] for r in pi)
plat_both=[r for r in plat_summary if r["PlatformPooledSignsConcordant"]!="NOT_BOTH_AVAILABLE"]
plat_conc=sum(r["PlatformPooledSignsConcordant"]=="YES" for r in plat_both)
small_counts=Counter(r["Classification"] for r in small)

summary=[
"PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT",
"STAGE 4E COMPLETE — PRIMARY ROBUSTNESS / SENSITIVITY ANALYSES","",
"Primary Stage4D results altered/replaced: 0",
"GSE295540 accessed: NO",
"Post-hoc sign flipping: NO","",
"LEAVE-ONE-COHORT-OUT",
f"Paper×gene groups eligible for LOCO (k>=3): {len(loco_summary)}",
f"Groups replicated with CI>0 in every LOCO run: {loco_all_rep}",
f"Groups reversed with CI<0 in every LOCO run: {loco_all_rev}",
"",
"PLATFORM STRATIFICATION",
f"Paper×gene groups with both microarray and RNA-seq pooled estimates: {len(plat_both)}",
f"Groups with concordant pooled signs across both platforms: {plat_conc} / {len(plat_both)}",
"",
"SMALL-HEALTHY-GROUP EXCLUSION (N_Healthy < 10)",
f"Paper×gene syntheses affected: {len(small)}",
]
for k in ["REPLICATED_CI_EXCLUDES_ZERO","POSITIVE_INCONCLUSIVE","REVERSED_CI_EXCLUDES_ZERO",
          "NEGATIVE_INCONCLUSIVE","NULL_OR_TIE"]:
    summary.append(f"  {k}: {small_counts[k]}")
summary += [
"",
"MODIFIED HARTUNG-KNAPP",
]
for k in ["REPLICATED_CI_EXCLUDES_ZERO","POSITIVE_INCONCLUSIVE","REVERSED_CI_EXCLUDES_ZERO",
          "NEGATIVE_INCONCLUSIVE","NULL_OR_TIE"]:
    summary.append(f"  {k}: {mhk_counts[k]}")
summary += [
"",
"PREDICTION INTERVALS",
]
for k in ["PI_POSITIVE","PI_CROSSES_ZERO","PI_NEGATIVE","UNAVAILABLE_K_LT_3"]:
    summary.append(f"  {k}: {pi_counts[k]}")
summary += [
"",
"AUC ATTENUATION AFTER EXCLUDING N_Healthy < 10",
f"Rows retained: {len(att_kept)} / 202",
f"Median published-minus-independent attenuation: {statistics.median(atts):.6f}",
f"Mean published-minus-independent attenuation: {statistics.mean(atts):.6f}",
"",
"NEXT GATE",
"Use Stage4D as the primary result and Stage4E as robustness evidence.",
"If these sensitivities remain coherent, Stage4F can generate manuscript-ready primary tables",
"and figures while GSE295540 remains sealed for the final temporal holdout."
]
(OUT/"STAGE4E_SUMMARY.txt").write_text("\n".join(summary)+"\n",encoding="utf-8")

manifest=[]
for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {"STAGE4E_FILE_MANIFEST.csv","Stage4E_RobustnessSensitivity_Bundle.zip","STAGE4E_COMPLETE.flag"}:
        manifest.append({"RelativePath":str(p.relative_to(OUT)),"Bytes":p.stat().st_size,"SHA256":sha(p)})
wcsv(OUT/"STAGE4E_FILE_MANIFEST.csv",manifest)

bundle=OUT/"Stage4E_RobustnessSensitivity_Bundle.zip"
if bundle.exists(): bundle.unlink()
with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p!=bundle: z.write(p,p.relative_to(OUT))
COMPLETE.write_text("Completed. Primary Stage4D unchanged. Sensitivity analyses complete. GSE295540 not accessed.\n",encoding="utf-8")

print("\n==========================================================")
print("STAGE 4E COMPLETE")
print("==========================================================")
print("\n".join(summary))
print("\nBundle:")
print(bundle)


'@ | Set-Content -Encoding UTF8 -Path $Py
$PyExe=$null
if(Get-Command py.exe -ErrorAction SilentlyContinue){$PyExe='py.exe'}
elseif(Get-Command python.exe -ErrorAction SilentlyContinue){$PyExe='python.exe'}
elseif(Get-Command python -ErrorAction SilentlyContinue){$PyExe='python'}
else{throw 'Python 3 was not found.'}
if($PyExe -eq 'py.exe'){
    & $PyExe -3 -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 4E syntax check failed.'}
    & $PyExe -3 $Py
}else{
    & $PyExe -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 4E syntax check failed.'}
    & $PyExe $Py
}
if($LASTEXITCODE -ne 0){throw "Stage 4E failed with exit code $LASTEXITCODE"}
Stop-Transcript
