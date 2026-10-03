#requires -Version 7.0
$ErrorActionPreference='Stop'
Set-StrictMode -Version Latest

$Root=$env:PSORIASIS_CORRECTION_ROOT
if (-not $Root) { throw 'Set PSORIASIS_CORRECTION_ROOT to a separate correction workspace.' }
$ScriptDir=Join-Path $Root '04_scripts\stage5'
$LogDir=Join-Path $Root 'logs'
New-Item -ItemType Directory -Force -Path $ScriptDir,$LogDir | Out-Null

$Log=Join-Path $LogDir ("stage5c_postholdout_updated_synthesis_"+(Get-Date -Format 'yyyyMMdd_HHmmss')+".log")
Start-Transcript -Path $Log

Write-Host '=========================================================='
Write-Host 'STAGE 5C: SECONDARY POST-HOLDOUT UPDATED SYNTHESIS'
Write-Host '=========================================================='
Write-Host 'Stage4D remains the primary pre-holdout synthesis.'
Write-Host 'Stage5B remains the standalone final temporal holdout.'
Write-Host 'This stage appends GSE295540 only in a separately labeled secondary synthesis.'
Write-Host ''

$Py=Join-Path $ScriptDir 'stage5c_postholdout_updated_synthesis.py'
@'

# CORRECTION PATCH 2026-09-07. Preserve original frozen results.

import csv, hashlib, math, statistics, zipfile
from collections import Counter, defaultdict
from pathlib import Path

import os
ROOT = Path(os.environ["PSORIASIS_CORRECTION_ROOT"])
S4D0=ROOT/"05_results"/"primary"/"stage4d0_synthesis_independence_freeze"
S4D=ROOT/"05_results"/"primary"/"stage4d_auc_attenuation_and_random_effects"
S4E=ROOT/"05_results"/"sensitivity"/"stage4e_primary_robustness"
S5B=ROOT/"05_results"/"primary"/"stage5b_final_temporal_holdout"
OUT=ROOT/"05_results"/"secondary"/"stage5c_postholdout_updated_synthesis"

PREINPUT=S4D0/"single_gene_synthesis_ready_metrics_FINAL_FROZEN_STAGE4D0.csv"
PREMETA=S4D/"paper_gene_random_effects_hedges_g_STAGE4D.csv"
PREMHK=S4E/"modified_hartung_knapp_sensitivity_STAGE4E.csv"
HOLD=S5B/"GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv"

for p in [PREINPUT,PREMETA,PREMHK,HOLD,
          S4D0/"STAGE4D0_COMPLETE.flag",S4D/"STAGE4D_COMPLETE.flag",
          S4E/"STAGE4E_COMPLETE.flag",S5B/"STAGE5B_COMPLETE.flag"]:
    if not p.exists() or p.stat().st_size==0:
        raise RuntimeError(f"Required input missing/empty: {p}")

OUT.mkdir(parents=True,exist_ok=True)
COMPLETE=OUT/"STAGE5C_COMPLETE.flag"
if COMPLETE.exists():
    raise RuntimeError(f"Stage5C already complete: {COMPLETE}")

def clean(x): return (x or "").strip().strip('"').strip("'")
def rcsv(p):
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
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
STAGE 5C — SECONDARY POST-HOLDOUT UPDATED SYNTHESIS

Status
Secondary analysis performed only after the standalone GSE295540 temporal-holdout result file
was completely written and frozen in Stage5B.

Primary analysis remains Stage4D and is never overwritten or relabeled.

Input
- Stage4D0 synthesis-independent pre-holdout cohort metrics
- Stage5B standalone GSE295540 holdout metrics

Update rule
For each of the 39 frozen PMID × GeneSymbol groups:
1. retain every Stage4D0 synthesis input unchanged;
2. append exactly one GSE295540 aligned Hedges-g estimate from Stage5B;
3. use the same frozen Stage2G orientation;
4. rerun the same REML random-effects Hedges-g model used in Stage4D;
5. calculate normal-theory 95% CI, tau^2, Q, I^2, and prediction interval;
6. calculate a modified Hartung-Knapp interval using the same Stage4E rule;
7. compare the resulting classification with the pre-holdout Stage4D classification.

No marker selection, cohort exclusion, post-hoc sign flip, or comparator remapping is allowed.
No grand meta-analysis across the 39 markers is performed because marker syntheses reuse cohorts.

This analysis is explicitly POST-HOLDOUT SECONDARY and cannot replace the standalone Stage5B
holdout result or the primary pre-holdout Stage4D synthesis.
"""
(OUT/"STAGE5C_POSTHOLDOUT_SYNTHESIS_PROTOCOL_FROZEN.txt").write_text(PROTOCOL,encoding="utf-8")

pre=rcsv(PREINPUT)
premeta=rcsv(PREMETA)
premhk=rcsv(PREMHK)
hold=rcsv(HOLD)

if len(pre)!=202: raise RuntimeError(f"Expected 202 Stage4D0 rows; got {len(pre)}")
if len(premeta)!=39: raise RuntimeError(f"Expected 39 Stage4D meta rows; got {len(premeta)}")
if len(premhk)!=39: raise RuntimeError(f"Expected 39 Stage4E mHK rows; got {len(premhk)}")
if len(hold)!=39: raise RuntimeError(f"Expected 39 Stage5B rows; got {len(hold)}")
if any(clean(r["Availability"])!="AVAILABLE" for r in hold):
    raise RuntimeError("Stage5C requires all 39 holdout marker rows available.")

def key(r): return (clean(r["PMID"]),clean(r["GeneSymbol"]))
pmeta={key(r):r for r in premeta}
pmhk={key(r):r for r in premhk}
hmap={key(r):r for r in hold}
if len(pmeta)!=39 or len(pmhk)!=39 or len(hmap)!=39:
    raise RuntimeError("Duplicate/missing PMID×GeneSymbol keys in pre/holdout tables.")

groups=defaultdict(list)
for r in pre:
    groups[key(r)].append(r)
if len(groups)!=39:
    raise RuntimeError(f"Expected 39 pre-holdout groups; found {len(groups)}")
if set(groups)!=set(hmap):
    raise RuntimeError("Stage4D0 and Stage5B PMID×GeneSymbol key sets differ.")

def reml_components(y,v,tau2):
    w=[1.0/(vi+tau2) for vi in v]
    sw=sum(w)
    mu=sum(wi*yi for wi,yi in zip(w,y))/sw
    q=sum(wi*(yi-mu)**2 for wi,yi in zip(w,y))
    obj=0.5*(sum(math.log(vi+tau2) for vi in v)+math.log(sw)+q)
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
    f=lambda t: reml_components(y,v,t)[0]
    f0=f(0.0)
    vy=statistics.variance(y) if len(y)>=2 else 0.0
    upper=max(1e-8,vy,statistics.mean(v),0.01)
    prev=f(upper)
    for _ in range(30):
        nxt=upper*4.0
        fn=f(nxt)
        if fn>=prev:
            upper = nxt  # correction: retain the rising endpoint in the bracket
            break
        upper=nxt; prev=fn
    th,obj=golden(f,0.0,upper)
    return 0.0 if f0<=obj+1e-10 else max(0.0,th)

def meta(y,se):
    k=len(y); v=[s*s for s in se]
    if k<2: raise RuntimeError("Stage5C groups must have >=2 estimates")
    if any((not math.isfinite(x) or x<=0) for x in v):
        raise RuntimeError(f"Invalid variance: {v}")
    wf=[1/x for x in v]
    swf=sum(wf)
    muf=sum(wi*yi for wi,yi in zip(wf,y))/swf
    Q=sum(wi*(yi-muf)**2 for wi,yi in zip(wf,y))
    df=k-1
    I2=max(0.0,(Q-df)/Q*100.0) if Q>0 else 0.0
    t2=tau_reml(y,v)
    _,mu,w,qre=reml_components(y,v,t2)
    sem=math.sqrt(1.0/sum(w))
    lo=mu-1.96*sem; hi=mu+1.96*sem
    pse=math.sqrt(t2+sem*sem)
    plo=mu-1.96*pse; phi=mu+1.96*pse
    return {"k":k,"mu":mu,"se":sem,"lo":lo,"hi":hi,"tau2":t2,"Q":Q,"df":df,"I2":I2,
            "plo":plo,"phi":phi,"weights":w,"qre":qre}

def cls(mu,lo,hi):
    if mu>0 and lo>0: return "REPLICATED_CI_EXCLUDES_ZERO"
    if mu>0: return "POSITIVE_INCONCLUSIVE"
    if mu<0 and hi<0: return "REVERSED_CI_EXCLUDES_ZERO"
    if mu<0: return "NEGATIVE_INCONCLUSIVE"
    return "NULL_OR_TIE"

def picls(lo,hi):
    if lo>0: return "PI_POSITIVE"
    if hi<0: return "PI_NEGATIVE"
    return "PI_CROSSES_ZERO"

# Exact 0.975 t criticals for df 1..30; only df <=7 expected here.
T975={1:12.7062047364,2:4.30265272991,3:3.18244630528,4:2.7764451052,
      5:2.57058183564,6:2.44691185114,7:2.36462425101,8:2.30600413503,
      9:2.26215716285,10:2.22813885196}

updated=[]
comparison=[]
for k,rs in sorted(groups.items()):
    hr=hmap[k]
    y=[float(r["HedgesG_Aligned"]) for r in rs] + [float(hr["HedgesG_Aligned"])]
    se=[float(r["HedgesG_SE"]) for r in rs] + [float(hr["HedgesG_SE"])]
    m=meta(y,se)
    primary_cls=cls(m["mu"],m["lo"],m["hi"])
    pi_class=picls(m["plo"],m["phi"])

    df=m["k"]-1
    if df not in T975:
        raise RuntimeError(f"No frozen t critical for df={df}")
    qstar=max(1.0,m["qre"]/df)
    se_mhk=math.sqrt(qstar/sum(m["weights"]))
    crit=T975[df]
    mhk_lo=m["mu"]-crit*se_mhk
    mhk_hi=m["mu"]+crit*se_mhk
    mhk_cls=cls(m["mu"],mhk_lo,mhk_hi)

    pm=pmeta[k]; ph=pmhk[k]
    pre_mu=float(pm["PooledHedgesG_REML"])
    pre_lo=float(pm["PooledHedgesG_CI95_Lower"])
    pre_hi=float(pm["PooledHedgesG_CI95_Upper"])
    pre_class=clean(pm["SynthesisClassification"])
    pre_mhk_class=clean(ph["mHK_Classification"])

    updated.append({
        "PMID":k[0],"GeneSymbol":k[1],
        "K_PreHoldout":len(rs),"K_PostHoldout":m["k"],
        "AppendedCohort":"GSE295540",
        "HoldoutHedgesG":hr["HedgesG_Aligned"],
        "HoldoutHedgesG_SE":hr["HedgesG_SE"],
        "UpdatedPooledHedgesG_REML":f"{m['mu']:.12g}",
        "UpdatedPooledHedgesG_SE":f"{m['se']:.12g}",
        "UpdatedCI95_Lower":f"{m['lo']:.12g}",
        "UpdatedCI95_Upper":f"{m['hi']:.12g}",
        "UpdatedTau2_REML":f"{m['tau2']:.12g}",
        "UpdatedCochranQ":f"{m['Q']:.12g}",
        "UpdatedQ_df":m["df"],
        "UpdatedI2_Percent":f"{m['I2']:.12g}",
        "UpdatedPredictionInterval_Lower":f"{m['plo']:.12g}",
        "UpdatedPredictionInterval_Upper":f"{m['phi']:.12g}",
        "UpdatedPredictionIntervalClassification":pi_class,
        "UpdatedPrimaryClassification":primary_cls,
        "Updated_mHK_qstar":f"{qstar:.12g}",
        "Updated_mHK_SE":f"{se_mhk:.12g}",
        "Updated_mHK_CI95_Lower":f"{mhk_lo:.12g}",
        "Updated_mHK_CI95_Upper":f"{mhk_hi:.12g}",
        "Updated_mHK_Classification":mhk_cls,
        "AnalysisLabel":"SECONDARY_POSTHOLDOUT_UPDATED_SYNTHESIS",
    })

    comparison.append({
        "PMID":k[0],"GeneSymbol":k[1],
        "PreHoldoutPooledG":f"{pre_mu:.12g}",
        "PostHoldoutPooledG":f"{m['mu']:.12g}",
        "ChangeInPooledG":f"{m['mu']-pre_mu:.12g}",
        "PreHoldoutPrimaryClassification":pre_class,
        "PostHoldoutPrimaryClassification":primary_cls,
        "PrimaryClassificationChanged":"YES" if primary_cls!=pre_class else "NO",
        "PreHoldout_mHK_Classification":pre_mhk_class,
        "PostHoldout_mHK_Classification":mhk_cls,
        "mHK_ClassificationChanged":"YES" if mhk_cls!=pre_mhk_class else "NO",
        "HoldoutSignMatchesPreHoldoutPooledG":"YES" if (float(hr["HedgesG_Aligned"])>0)==(pre_mu>0) else "NO",
        "PreHoldoutI2_Percent":pm["I2_Percent"],
        "PostHoldoutI2_Percent":f"{m['I2']:.12g}",
        "I2_Change":f"{m['I2']-float(pm['I2_Percent']):.12g}",
    })

wcsv(OUT/"secondary_postholdout_random_effects_STAGE5C.csv",updated)
wcsv(OUT/"pre_vs_post_holdout_synthesis_comparison_STAGE5C.csv",comparison)

primary_counts=Counter(r["UpdatedPrimaryClassification"] for r in updated)
mhk_counts=Counter(r["Updated_mHK_Classification"] for r in updated)
pi_counts=Counter(r["UpdatedPredictionIntervalClassification"] for r in updated)
prim_changes=[r for r in comparison if r["PrimaryClassificationChanged"]=="YES"]
mhk_changes=[r for r in comparison if r["mHK_ClassificationChanged"]=="YES"]
all_sign=sum(r["HoldoutSignMatchesPreHoldoutPooledG"]=="YES" for r in comparison)
pooled_changes=[float(r["ChangeInPooledG"]) for r in comparison]
i2_changes=[float(r["I2_Change"]) for r in comparison]

lowest=sorted(updated,key=lambda r:float(r["UpdatedPooledHedgesG_REML"]))

summary=[
"PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT",
"STAGE 5C COMPLETE — SECONDARY POST-HOLDOUT UPDATED SYNTHESIS","",
"Primary Stage4D analysis overwritten/replaced: NO",
"Standalone Stage5B holdout overwritten/replaced: NO",
"Post-hoc sign flipping: NO",
"Performance-driven exclusions: NO",
"Grand meta-analysis across markers: NO","",
f"Updated paper×gene syntheses: {len(updated)} / 39",
f"GSE295540 effect sign matched pre-holdout pooled-g sign: {all_sign} / 39",
f"Primary classification changes after adding holdout: {len(prim_changes)} / 39",
f"mHK classification changes after adding holdout: {len(mhk_changes)} / 39",
f"Median change in pooled Hedges g: {statistics.median(pooled_changes):.6f}",
f"Median change in I2: {statistics.median(i2_changes):.2f} percentage points","",
"UPDATED NORMAL-THEORY REML CLASSIFICATIONS",
]
for x in ["REPLICATED_CI_EXCLUDES_ZERO","POSITIVE_INCONCLUSIVE","REVERSED_CI_EXCLUDES_ZERO",
          "NEGATIVE_INCONCLUSIVE","NULL_OR_TIE"]:
    summary.append(f"  {x}: {primary_counts[x]}")
summary += ["","UPDATED MODIFIED HARTUNG-KNAPP CLASSIFICATIONS"]
for x in ["REPLICATED_CI_EXCLUDES_ZERO","POSITIVE_INCONCLUSIVE","REVERSED_CI_EXCLUDES_ZERO",
          "NEGATIVE_INCONCLUSIVE","NULL_OR_TIE"]:
    summary.append(f"  {x}: {mhk_counts[x]}")
summary += ["","UPDATED PREDICTION INTERVALS"]
for x in ["PI_POSITIVE","PI_CROSSES_ZERO","PI_NEGATIVE"]:
    summary.append(f"  {x}: {pi_counts[x]}")

summary += ["","LOWEST UPDATED POOLED EFFECTS"]
for r in lowest[:8]:
    summary.append(
        f"  PMID {r['PMID']} / {r['GeneSymbol']}: "
        f"g={float(r['UpdatedPooledHedgesG_REML']):.6f} "
        f"[{float(r['UpdatedCI95_Lower']):.6f}, {float(r['UpdatedCI95_Upper']):.6f}], "
        f"mHK={r['Updated_mHK_Classification']}, PI={r['UpdatedPredictionIntervalClassification']}"
    )

summary += ["","CLASSIFICATION CHANGES"]
if prim_changes:
    for r in prim_changes:
        summary.append(
            f"  PRIMARY PMID {r['PMID']} / {r['GeneSymbol']}: "
            f"{r['PreHoldoutPrimaryClassification']} -> {r['PostHoldoutPrimaryClassification']}"
        )
else:
    summary.append("  Primary REML: NONE")
if mhk_changes:
    for r in mhk_changes:
        summary.append(
            f"  mHK PMID {r['PMID']} / {r['GeneSymbol']}: "
            f"{r['PreHoldout_mHK_Classification']} -> {r['PostHoldout_mHK_Classification']}"
        )
else:
    summary.append("  Modified Hartung-Knapp: NONE")

summary += [
    "",
    "INTERPRETIVE STATUS",
    "Stage4D remains the primary pre-holdout synthesis.",
    "Stage5B remains the standalone final temporal holdout.",
    "Stage5C is secondary evidence describing how the pooled estimates change after the holdout is appended.",
    "",
    "NEXT GATE",
    "Stage5D can generate final manuscript-ready integrated tables/figures that preserve",
    "the three-layer evidence structure: primary pre-holdout, standalone temporal holdout,",
    "and secondary post-holdout updated synthesis."
]
(OUT/"STAGE5C_SUMMARY.txt").write_text("\n".join(summary)+"\n",encoding="utf-8")

manifest=[]
for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {"STAGE5C_FILE_MANIFEST.csv","Stage5C_PostHoldoutUpdatedSynthesis_Bundle.zip","STAGE5C_COMPLETE.flag"}:
        manifest.append({"RelativePath":str(p.relative_to(OUT)),"Bytes":p.stat().st_size,"SHA256":sha(p)})
wcsv(OUT/"STAGE5C_FILE_MANIFEST.csv",manifest)

bundle=OUT/"Stage5C_PostHoldoutUpdatedSynthesis_Bundle.zip"
if bundle.exists(): bundle.unlink()
with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p!=bundle:
            z.write(p,p.relative_to(OUT))

COMPLETE.write_text(
    "Completed. Secondary post-holdout updated synthesis only. Stage4D and Stage5B unchanged.\n",
    encoding="utf-8"
)

print("\n==========================================================")
print("STAGE 5C COMPLETE")
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
    if($LASTEXITCODE -ne 0){throw 'Stage 5C syntax check failed.'}
    & $PyExe -3 $Py
}else{
    & $PyExe -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 5C syntax check failed.'}
    & $PyExe $Py
}
if($LASTEXITCODE -ne 0){throw "Stage 5C failed with exit code $LASTEXITCODE"}

Stop-Transcript
