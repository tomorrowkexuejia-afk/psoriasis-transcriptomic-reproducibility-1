#requires -Version 7.0
$ErrorActionPreference='Stop'
Set-StrictMode -Version Latest

$Root=$env:PSORIASIS_CORRECTION_ROOT
if (-not $Root) { throw 'Set PSORIASIS_CORRECTION_ROOT to a separate correction workspace.' }
$ScriptDir=Join-Path $Root '04_scripts\stage5'
$LogDir=Join-Path $Root 'logs'
New-Item -ItemType Directory -Force -Path $ScriptDir,$LogDir | Out-Null

$Log=Join-Path $LogDir ("stage5b_final_temporal_holdout_"+(Get-Date -Format 'yyyyMMdd_HHmmss')+".log")
Start-Transcript -Path $Log

Write-Host '=========================================================='
Write-Host 'STAGE 5B: FINAL TEMPORAL HOLDOUT — GSE295540'
Write-Host '=========================================================='
Write-Host 'This run unseals the frozen GSE295540 count matrix under Stage5A.'
Write-Host 'The matrix is streamed once; Stage4D is NOT updated.'
Write-Host 'No post-hoc sign flipping or performance-driven exclusions.'
Write-Host ''

$Py=Join-Path $ScriptDir 'stage5b_final_temporal_holdout.py'
@'

# CORRECTION PATCH 2026-09-07. Preserve original frozen results.

import csv, gzip, hashlib, math, re, statistics, zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import os
ROOT = Path(os.environ["PSORIASIS_CORRECTION_ROOT"])
S5A=ROOT/"00_protocol"/"stage5a_final_temporal_holdout_freeze"
RAW=ROOT/"02_data"/"raw_geo"/"GSE295540"
OUT=ROOT/"05_results"/"primary"/"stage5b_final_temporal_holdout"

COUNTS=RAW/"GSE295540_raw_counts_All_samples.csv.gz"
SAMPLEMAP=S5A/"GSE295540_sample_phenotype_map_FINAL_FROZEN.csv"
TARGETS=S5A/"GSE295540_paper_gene_targets_FINAL_FROZEN.csv"
NOM=S5A/"GSE295540_target_nomenclature_FINAL_FROZEN.csv"
AUDIT5A=S5A/"STAGE5A_HOLDOUT_FREEZE_AUDIT.csv"
PROTOCOL=S5A/"STAGE5A_GSE295540_FINAL_TEMPORAL_HOLDOUT_PROTOCOL_FROZEN.txt"
COMPLETE5A=S5A/"STAGE5A_COMPLETE.flag"

for p in [COUNTS,SAMPLEMAP,TARGETS,NOM,AUDIT5A,PROTOCOL,COMPLETE5A]:
    if not p.exists() or p.stat().st_size==0:
        raise RuntimeError(f"Required frozen Stage5A input missing/empty: {p}")

OUT.mkdir(parents=True,exist_ok=True)
COMPLETE=OUT/"STAGE5B_COMPLETE.flag"
UNSEAL=OUT/"STAGE5B_UNSEAL_STARTED.flag"
if COMPLETE.exists():
    raise RuntimeError(f"Stage5B already complete/frozen: {COMPLETE}")

def clean(x): return (x or "").strip().strip('"').strip("'")
def norm(x): return re.sub(r"[^a-z0-9]+","",clean(x).lower())
def rcsv(p):
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
def wcsv(p,rows,fields=None):
    rows=list(rows)
    if fields is None: fields=list(rows[0].keys()) if rows else []
    with p.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        if rows: w.writerows(rows)
def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def sgn(x,eps=1e-15):
    return 1 if x>eps else -1 if x<-eps else 0

# ---------------- Frozen identity checks before unseal ----------------
a5=rcsv(AUDIT5A)
frozen_sha=None
for r in a5:
    if clean(r.get("Item"))=="Compressed SHA256":
        frozen_sha=clean(r.get("Value"))
    if clean(r.get("Item"))=="Raw-count SHA256":
        frozen_sha=clean(r.get("Value"))
if not frozen_sha:
    raise RuntimeError("Could not find frozen GSE295540 compressed SHA256 in Stage5A audit.")

current_sha=sha256(COUNTS)
if current_sha.lower()!=frozen_sha.lower():
    raise RuntimeError(f"HOLDOUT FILE HASH MISMATCH. Frozen={frozen_sha} Current={current_sha}")

samplemap=rcsv(SAMPLEMAP)
targets=rcsv(TARGETS)
nom=rcsv(NOM)
if len(samplemap)!=13: raise RuntimeError(f"Frozen sample map has {len(samplemap)} rows; expected 13.")
if len(targets)!=39: raise RuntimeError(f"Frozen target table has {len(targets)} rows; expected 39.")
if len({(r["PMID"],r["GeneSymbol"]) for r in targets})!=39:
    raise RuntimeError("Frozen target table has duplicate PMID×GeneSymbol keys.")

# Frozen sample alias mapping.
alias_to_gsm={}
sample_by_gsm={}
for r in samplemap:
    gsm=clean(r["GSM"])
    sample_by_gsm[gsm]=r
    aliases=[clean(x) for x in clean(r["FrozenExpressionColumnAliases"]).split("|") if clean(x)]
    for a in aliases:
        n=norm(a)
        if not n: continue
        if n in alias_to_gsm and alias_to_gsm[n]!=gsm:
            raise RuntimeError(f"Frozen sample alias ambiguity: {a!r} maps to {alias_to_gsm[n]} and {gsm}")
        alias_to_gsm[n]=gsm

# Frozen target nomenclature maps.
symbol_to_gid=defaultdict(set)
ens_to_gid=defaultdict(set)
gid_to_nom={}
for r in nom:
    gid=clean(r["NCBIGeneID"])
    gid_to_nom[gid]=r
    for a in clean(r["FrozenSymbolAliases"]).split("|"):
        a=clean(a)
        if a: symbol_to_gid[a.upper()].add(gid)
    for e in clean(r["FrozenEnsemblGeneIDs"]).split("|"):
        e=clean(e)
        if e: ens_to_gid[e.upper()].add(gid)

gid_to_current={clean(r["NCBIGeneID"]):clean(r["NCBICurrentSymbol"]) for r in nom}

# Only aliases unique among the frozen target set are usable.
# Correction release: deposited matrix has all targets as current symbols.
# Do not aggregate aliases such as the separate MIF gene into S100A9.
symbol_unique={clean(r["NCBICurrentSymbol"]).upper():clean(r["NCBIGeneID"]) for r in nom}
ens_unique={a:next(iter(gs)) for a,gs in ens_to_gid.items() if len(gs)==1}
target_gids={clean(r["NCBIGeneID"]) for r in targets}

# ---------------- Statistical functions frozen to Stage4A definitions ----------------
def midrank(x):
    n=len(x)
    order=sorted(range(n),key=lambda i:x[i])
    z=[x[i] for i in order]
    ranks=[0.0]*n
    i=0
    while i<n:
        j=i
        while j<n and z[j]==z[i]: j+=1
        mr=0.5*(i+j-1)+1
        for k in range(i,j): ranks[k]=mr
        i=j
    out=[0.0]*n
    for pos,orig in enumerate(order): out[orig]=ranks[pos]
    return out

def sample_var(a):
    if len(a)<2: return 0.0
    mu=statistics.mean(a)
    return sum((x-mu)**2 for x in a)/(len(a)-1)

def auc_delong(cases,ctrls):
    m=len(cases); n=len(ctrls)
    tx=midrank(cases); ty=midrank(ctrls); tz=midrank(cases+ctrls)
    auc=(sum(tz[:m])/m-(m+1)/2)/n
    v01=[(tz[i]-tx[i])/n for i in range(m)]
    v10=[1-(tz[m+j]-ty[j])/m for j in range(n)]
    var=max(0.0,sample_var(v01)/m+sample_var(v10)/n)
    se=math.sqrt(var)
    lo=max(0.0,auc-1.96*se); hi=min(1.0,auc+1.96*se)
    # Independent pairwise implementation check.
    auc2=sum(1 if p>h else .5 if p==h else 0 for p in cases for h in ctrls)/(m*n)
    if abs(auc-auc2)>1e-10:
        raise RuntimeError(f"AUC self-check failed: DeLong rank={auc} pairwise={auc2}")
    return auc,se,lo,hi

def hedges(cases,ctrls):
    n1=len(cases); n0=len(ctrls)
    m1=statistics.mean(cases); m0=statistics.mean(ctrls)
    s1=statistics.stdev(cases); s0=statistics.stdev(ctrls)
    df=n1+n0-2
    pv=((n1-1)*s1*s1+(n0-1)*s0*s0)/df
    if pv==0:
        if m1==m0:
            return 0.0,0.0,0.0,0.0,m1,m0
        raise RuntimeError("Zero pooled SD with nonzero mean difference in holdout.")
    d=(m1-m0)/math.sqrt(pv)
    J=1-3/(4*(n1+n0)-9)
    g=J*d
    vg=J*J*((n1+n0)/(n1*n0)+d*d/(2*df))
    se=math.sqrt(max(0.0,vg))
    return g,se,g-1.96*se,g+1.96*se,m1,m0

# ---------------- Unseal state ----------------
recovery=UNSEAL.exists()
if not recovery:
    UNSEAL.write_text(
        f"Stage5B holdout unseal began {datetime.now(timezone.utc).isoformat()}\n"
        f"Frozen compressed SHA256: {frozen_sha}\n"
        "Expression access occurs after this marker.\n",
        encoding="utf-8"
    )
else:
    # Do not pretend the holdout is still sealed. A deterministic recovery rerun is allowed.
    with UNSEAL.open("a",encoding="utf-8") as f:
        f.write(f"Deterministic recovery rerun under unchanged frozen protocol: {datetime.now(timezone.utc).isoformat()}\n")

# ---------------- SINGLE STREAMING PASS THROUGH HOLDOUT COUNTS ----------------
# Buffer first rows only for identifier-column selection; they are processed before continuing.
BUFFER_N=1000
with gzip.open(COUNTS,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
    reader=csv.reader(f)
    try:
        header=next(reader)
    except StopIteration:
        raise RuntimeError("GSE295540 count matrix is empty.")
    header=[clean(x) for x in header]
    if len(header)<14:
        raise RuntimeError(f"Count matrix has only {len(header)} columns; expected >=14.")

    # Resolve exactly 13 sample columns by frozen aliases.
    sample_col={}
    col_sample={}
    for i,h in enumerate(header):
        n=norm(h)
        if n in alias_to_gsm:
            gsm=alias_to_gsm[n]
            if gsm in sample_col:
                raise RuntimeError(f"Multiple count columns map to frozen sample {gsm}: {header[sample_col[gsm]]!r} and {h!r}")
            sample_col[gsm]=i
            col_sample[i]=gsm

    missing=sorted(set(sample_by_gsm)-set(sample_col))
    if missing:
        raise RuntimeError(
            "SAMPLE COLUMN RESOLUTION ABORT before phenotype-dependent metrics. "
            f"Matched {len(sample_col)}/13; missing={missing}; header={header}"
        )
    if len(sample_col)!=13:
        raise RuntimeError(f"Expected exactly 13 resolved sample columns; got {len(sample_col)}")

    non_sample_idx=[i for i in range(len(header)) if i not in col_sample]
    if not non_sample_idx:
        raise RuntimeError("No non-sample identifier/annotation columns remain.")

    buffer=[]
    for _ in range(BUFFER_N):
        try: row=next(reader)
        except StopIteration: break
        if len(row)!=len(header):
            raise RuntimeError(f"Malformed count row in identifier buffer: {len(row)} values vs {len(header)} header columns.")
        buffer.append(row)
    if not buffer:
        raise RuntimeError("Count matrix contains header but no data rows.")

    # Deterministically select identifier column using frozen target mappings only.
    ens_re=re.compile(r"^ENSG\d+(?:\.\d+)?$",re.I)
    candidates=[]
    for idx in non_sample_idx:
        vals=[clean(r[idx]) for r in buffer if clean(r[idx])]
        if not vals: continue
        ens_hits=set()
        sym_hits=set()
        ens_pattern=sum(bool(ens_re.match(v)) for v in vals)
        for v in vals:
            base=re.sub(r"\.\d+$","",v).upper()
            if base in ens_unique: ens_hits.add(ens_unique[base])
            if v.upper() in symbol_unique: sym_hits.add(symbol_unique[v.upper()])
        modes=[]
        if ens_hits: modes.append(("ENSEMBL",len(ens_hits)))
        if sym_hits: modes.append(("SYMBOL",len(sym_hits)))
        for mode,hits in modes:
            hn=norm(header[idx])
            priority=0
            if mode=="ENSEMBL":
                if "ensembl" in hn: priority=4
                elif hn in {"geneid","geneidentifier","id"} or "geneid" in hn: priority=3
                elif "gene" in hn: priority=2
            else:
                if "symbol" in hn: priority=4
                elif "genename" in hn or hn=="gene": priority=3
                elif "name" in hn: priority=2
            candidates.append({
                "idx":idx,"header":header[idx],"mode":mode,
                "target_hits":hits,"header_priority":priority,
                "ens_pattern_rate":ens_pattern/len(vals)
            })

    if not candidates:
        diag=[header[i] for i in non_sample_idx]
        raise RuntimeError(
            "IDENTIFIER RESOLUTION ABORT before phenotype-dependent metrics. "
            f"No identifier column matched any frozen target symbol/Ensembl mapping. Non-sample columns={diag}"
        )

    candidates.sort(key=lambda c:(c["target_hits"],c["header_priority"]),reverse=True)
    best=candidates[0]
    tied=[
        c for c in candidates
        if (c["target_hits"],c["header_priority"])==(best["target_hits"],best["header_priority"])
        and (c["idx"],c["mode"])!=(best["idx"],best["mode"])
    ]
    if tied:
        raise RuntimeError(
            "IDENTIFIER RESOLUTION ABORT: multiple equally preferred frozen-compatible identifier mappings. "
            f"Best={best}; ties={tied}"
        )
    id_idx=best["idx"]; id_mode=best["mode"]

    # State accumulated in the single stream.
    library={gsm:0.0 for gsm in sample_by_gsm}
    target_counts={gid:{gsm:0.0 for gsm in sample_by_gsm} for gid in target_gids}
    target_row_counts=Counter()
    total_rows=0
    nonblank_ids=0
    mixed_identifier_examples=[]

    def process(row):
        global total_rows,nonblank_ids
        total_rows+=1
        rid=clean(row[id_idx])
        if not rid:
            return
        nonblank_ids+=1

        # Counts are validated and summed for ALL rows before target filtering.
        vals={}
        for gsm,idx in sample_col.items():
            raw=clean(row[idx])
            if raw=="":
                raise RuntimeError(f"Blank count at data row {total_rows}, sample {gsm}")
            try: v=float(raw)
            except Exception:
                raise RuntimeError(f"Non-numeric count {raw!r} at data row {total_rows}, sample {gsm}")
            if not math.isfinite(v) or v<0:
                raise RuntimeError(f"Invalid count {v} at data row {total_rows}, sample {gsm}")
            vals[gsm]=v
            library[gsm]+=v

        gid=None
        if id_mode=="ENSEMBL":
            if not ens_re.match(rid):
                if len(mixed_identifier_examples)<10:
                    mixed_identifier_examples.append(rid)
            base=re.sub(r"\.\d+$","",rid).upper()
            gid=ens_unique.get(base)
        else:
            if ens_re.match(rid):
                if len(mixed_identifier_examples)<10:
                    mixed_identifier_examples.append(rid)
            gid=symbol_unique.get(rid.upper())

        if gid in target_counts:
            target_row_counts[gid]+=1
            for gsm,v in vals.items():
                target_counts[gid][gsm]+=v

    for row in buffer: process(row)
    for row in reader:
        if len(row)!=len(header):
            raise RuntimeError(f"Malformed count row {total_rows+1}: {len(row)} values vs {len(header)} header columns.")
        process(row)

# After full single-pass structural read, enforce identifier-mode consistency.
if mixed_identifier_examples:
    raise RuntimeError(
        "IDENTIFIER RESOLUTION ABORT: selected identifier column appears mixed/incompatible "
        f"with frozen {id_mode} mode. Examples={mixed_identifier_examples}"
    )
if total_rows==0 or nonblank_ids==0:
    raise RuntimeError("No data rows/nonblank identifiers parsed.")
if any(v<=0 for v in library.values()):
    raise RuntimeError(f"At least one holdout library size is <=0: {library}")

# Mapping audit before performance calculations.
mapaudit=[]
for gid in sorted(target_gids,key=lambda x:int(x) if x.isdigit() else x):
    mapaudit.append({
        "NCBIGeneID":gid,
        "NCBICurrentSymbol":gid_to_current.get(gid,""),
        "RawRowsMapped":target_row_counts[gid],
        "Availability":"AVAILABLE" if target_row_counts[gid]>0 else "UNAVAILABLE_IDENTIFIER_NOT_FOUND",
        "IdentifierMode":id_mode,
        "IdentifierColumn":header[id_idx],
    })
wcsv(OUT/"GSE295540_target_mapping_audit_STAGE5B.csv",mapaudit)

structure=[{
    "FrozenCompressedSHA256":frozen_sha,
    "CurrentCompressedSHA256":current_sha,
    "HashMatch":"YES",
    "UnsealRecoveryRerun":"YES" if recovery else "NO",
    "HeaderColumnCount":len(header),
    "ResolvedSampleColumns":len(sample_col),
    "IdentifierColumn":header[id_idx],
    "IdentifierColumnIndexZeroBased":id_idx,
    "IdentifierMode":id_mode,
    "IdentifierSelectionTargetHitsInFirst1000Rows":best["target_hits"],
    "IdentifierHeaderPriority":best["header_priority"],
    "TotalDepositedDataRows":total_rows,
    "TargetGenesAvailable":sum(target_row_counts[g]>0 for g in target_gids),
    "TargetGenesFrozen":len(target_gids),
}]
wcsv(OUT/"GSE295540_matrix_structure_audit_STAGE5B.csv",structure)

librows=[]
for gsm in sorted(library):
    sm=sample_by_gsm[gsm]
    librows.append({
        "GSM":gsm,
        "Title":sm["Title"],
        "Phenotype":sm["FinalPhenotype"],
        "FrozenMatchedHeader":header[sample_col[gsm]],
        "FullLibrarySize":f"{library[gsm]:.12g}",
    })
wcsv(OUT/"GSE295540_library_sizes_STAGE5B.csv",librows)

# ---------------- Frozen normalization ----------------
expr={}
for gid in target_gids:
    if target_row_counts[gid]==0: continue
    expr[gid]={}
    for gsm in sample_by_gsm:
        c=target_counts[gid][gsm]
        cpm=1_000_000.0*c/library[gsm]
        expr[gid][gsm]=math.log2(cpm+1.0)

# ---------------- Per paper×gene holdout metrics ----------------
results=[]
metric_cache={}
for t in targets:
    pmid=clean(t["PMID"]); gene=clean(t["GeneSymbol"]); gid=clean(t["NCBIGeneID"])
    ori=clean(t["FrozenOrientation"])
    cachekey=(gid,ori)

    if gid not in expr:
        results.append({
            "PMID":pmid,"GeneSymbol":gene,"NCBIGeneID":gid,"FrozenOrientation":ori,
            "Availability":"UNAVAILABLE_IDENTIFIER_NOT_FOUND",
            "N_Psoriasis":"7","N_Healthy":"6",
            "DirectionalAUC":"","AUC_SE_DeLong":"","AUC_CI95_Lower":"","AUC_CI95_Upper":"",
            "HedgesG_Aligned":"","HedgesG_SE":"","HedgesG_CI95_Lower":"","HedgesG_CI95_Upper":"",
            "OrientedPsoriasisMean":"","OrientedHealthyMean":"","OrientedMeanDifference":"",
            "DirectionConcordantWithPublished":"",
            "PublishedComparatorAUC":t["HoldoutPublishedComparatorAUC"],
            "PublishedMinusHoldoutAUCAttenuation":"",
            "PreHoldoutMedianIndependentDirectionalAUC":t["PreHoldoutMedianIndependentDirectionalAUC"],
            "HoldoutMinusPreHoldoutMedianIndependentAUC":"",
            "PreHoldoutPooledHedgesG":t["PreHoldoutPooledHedgesG"],
            "HoldoutVsPreHoldoutEmpiricalSignConcordance":"",
            "PostHocSignFlip":"NO",
        })
        continue

    if cachekey not in metric_cache:
        sign=1 if ori=="UP_IN_PSORIASIS" else -1
        cases=[]; ctrls=[]
        raw_case=[]; raw_ctrl=[]
        for gsm,sm in sample_by_gsm.items():
            v=expr[gid][gsm]
            ov=sign*v
            if sm["FinalPhenotype"]=="PSORIASIS_LESIONAL_SKIN":
                cases.append(ov); raw_case.append(v)
            elif sm["FinalPhenotype"]=="HEALTHY_CONTROL_SKIN":
                ctrls.append(ov); raw_ctrl.append(v)
            else:
                raise RuntimeError(f"Unexpected frozen phenotype for {gsm}: {sm['FinalPhenotype']}")
        if len(cases)!=7 or len(ctrls)!=6:
            raise RuntimeError(f"{gene}: expected 7/6 phenotype values, got {len(cases)}/{len(ctrls)}")
        auc,ase,alo,ahi=auc_delong(cases,ctrls)
        g,gse,glo,ghi,omc,omh=hedges(cases,ctrls)
        metric_cache[cachekey]={
            "auc":auc,"ase":ase,"alo":alo,"ahi":ahi,
            "g":g,"gse":gse,"glo":glo,"ghi":ghi,
            "omc":omc,"omh":omh,
            "dir":"YES" if omc>omh else "NO" if omc<omh else "TIE",
        }

    m=metric_cache[cachekey]
    pub=float(t["HoldoutPublishedComparatorAUC"])
    preauc=float(t["PreHoldoutMedianIndependentDirectionalAUC"])
    preg=float(t["PreHoldoutPooledHedgesG"])
    empirical="YES" if sgn(m["g"])==sgn(preg) else "TIE" if sgn(m["g"])==0 or sgn(preg)==0 else "NO"

    results.append({
        "PMID":pmid,"GeneSymbol":gene,"NCBIGeneID":gid,"FrozenOrientation":ori,
        "Availability":"AVAILABLE",
        "N_Psoriasis":"7","N_Healthy":"6",
        "DirectionalAUC":f"{m['auc']:.12g}",
        "AUC_SE_DeLong":f"{m['ase']:.12g}",
        "AUC_CI95_Lower":f"{m['alo']:.12g}",
        "AUC_CI95_Upper":f"{m['ahi']:.12g}",
        "HedgesG_Aligned":f"{m['g']:.12g}",
        "HedgesG_SE":f"{m['gse']:.12g}",
        "HedgesG_CI95_Lower":f"{m['glo']:.12g}",
        "HedgesG_CI95_Upper":f"{m['ghi']:.12g}",
        "OrientedPsoriasisMean":f"{m['omc']:.12g}",
        "OrientedHealthyMean":f"{m['omh']:.12g}",
        "OrientedMeanDifference":f"{m['omc']-m['omh']:.12g}",
        "DirectionConcordantWithPublished":m["dir"],
        "PublishedComparatorAUC":f"{pub:.12g}",
        "PublishedMinusHoldoutAUCAttenuation":f"{pub-m['auc']:.12g}",
        "PreHoldoutMedianIndependentDirectionalAUC":f"{preauc:.12g}",
        "HoldoutMinusPreHoldoutMedianIndependentAUC":f"{m['auc']-preauc:.12g}",
        "PreHoldoutPooledHedgesG":f"{preg:.12g}",
        "HoldoutVsPreHoldoutEmpiricalSignConcordance":empirical,
        "PostHocSignFlip":"NO",
    })

if len(results)!=39: raise RuntimeError(f"Expected 39 holdout result rows; got {len(results)}")
wcsv(OUT/"GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv",results)

# Deduplicated gene×orientation metrics.
dedup={}
for r in results:
    k=(r["GeneSymbol"],r["FrozenOrientation"])
    if k not in dedup:
        dedup[k]={
            "GeneSymbol":r["GeneSymbol"],
            "NCBIGeneID":r["NCBIGeneID"],
            "FrozenOrientation":r["FrozenOrientation"],
            "Availability":r["Availability"],
            "DirectionalAUC":r["DirectionalAUC"],
            "AUC_CI95_Lower":r["AUC_CI95_Lower"],
            "AUC_CI95_Upper":r["AUC_CI95_Upper"],
            "HedgesG_Aligned":r["HedgesG_Aligned"],
            "HedgesG_CI95_Lower":r["HedgesG_CI95_Lower"],
            "HedgesG_CI95_Upper":r["HedgesG_CI95_Upper"],
            "DirectionConcordantWithPublished":r["DirectionConcordantWithPublished"],
            "PMIDs":set(),
        }
    else:
        # Duplicates must have identical expression-derived metrics.
        for c in ["Availability","DirectionalAUC","AUC_CI95_Lower","AUC_CI95_Upper",
                  "HedgesG_Aligned","HedgesG_CI95_Lower","HedgesG_CI95_Upper",
                  "DirectionConcordantWithPublished"]:
            if dedup[k][c]!=r[c]:
                raise RuntimeError(f"Duplicate gene×orientation metric inconsistency for {k}, field {c}")
    dedup[k]["PMIDs"].add(r["PMID"])

dedup_rows=[]
for k,d in sorted(dedup.items()):
    d=dict(d)
    d["PMIDs"]="|".join(sorted(d["PMIDs"]))
    dedup_rows.append(d)
if len(dedup_rows)!=36:
    raise RuntimeError(f"Expected 36 unique gene×orientation rows; got {len(dedup_rows)}")
wcsv(OUT/"GSE295540_unique_gene_orientation_holdout_results_STAGE5B.csv",dedup_rows)

# ---------------- Summary ----------------
avail=[r for r in results if r["Availability"]=="AVAILABLE"]
unavail=[r for r in results if r["Availability"]!="AVAILABLE"]
unique_avail=[r for r in dedup_rows if r["Availability"]=="AVAILABLE"]

def auc_bucket(a):
    if a>0.5: return "GT_0.5"
    if a<0.5: return "LT_0.5"
    return "EQ_0.5"

row_buckets=Counter(auc_bucket(float(r["DirectionalAUC"])) for r in avail)
uniq_buckets=Counter(auc_bucket(float(r["DirectionalAUC"])) for r in unique_avail)
aucs=[float(r["DirectionalAUC"]) for r in avail]
gs=[float(r["HedgesG_Aligned"]) for r in avail]
atts=[float(r["PublishedMinusHoldoutAUCAttenuation"]) for r in avail]
deltas=[float(r["HoldoutMinusPreHoldoutMedianIndependentAUC"]) for r in avail]
emp=Counter(r["HoldoutVsPreHoldoutEmpiricalSignConcordance"] for r in avail)

tl=[r for r in results if r["GeneSymbol"]=="TLN1"]
tl_line="TLN1 unavailable"
if tl and tl[0]["Availability"]=="AVAILABLE":
    x=tl[0]
    tl_line=(
        f"TLN1 holdout: AUC={float(x['DirectionalAUC']):.6f}; "
        f"g={float(x['HedgesG_Aligned']):.6f}; "
        f"published-direction concordance={x['DirectionConcordantWithPublished']}; "
        f"pre-holdout empirical-sign concordance={x['HoldoutVsPreHoldoutEmpiricalSignConcordance']}"
    )

summary=[
"PSORIASIS TRANSCRIPTOMIC BIOMARKER REPRODUCIBILITY AUDIT",
"STAGE 5B COMPLETE — FINAL TEMPORAL HOLDOUT GSE295540","",
f"Holdout file SHA-256 matched Stage5A freeze: YES",
f"Recovery rerun after prior unseal marker: {'YES' if recovery else 'NO'}",
"Stage4D synthesis updated/refit: NO",
"Post-hoc sign flipping: NO",
"Performance-driven exclusions: NO","",
"MATRIX / MAPPING",
f"Resolved frozen sample columns: 13 / 13",
f"Identifier column: {header[id_idx]}",
f"Identifier mode: {id_mode}",
f"Deposited data rows streamed: {total_rows}",
f"Unique target genes available: {sum(target_row_counts[g]>0 for g in target_gids)} / {len(target_gids)}",
f"Paper×gene rows available: {len(avail)} / 39",
f"Paper×gene rows unavailable: {len(unavail)}","",
"PAPER×GENE HOLDOUT RESULTS",
f"Directional AUC >0.5: {row_buckets['GT_0.5']} / {len(avail)}",
f"Directional AUC =0.5: {row_buckets['EQ_0.5']} / {len(avail)}",
f"Directional AUC <0.5: {row_buckets['LT_0.5']} / {len(avail)}",
]
if aucs:
    summary += [
        f"Median directional AUC: {statistics.median(aucs):.6f}",
        f"AUC range: {min(aucs):.6f} to {max(aucs):.6f}",
        f"Median aligned Hedges g: {statistics.median(gs):.6f}",
        f"Median published-minus-holdout AUC attenuation: {statistics.median(atts):.6f}",
        f"Mean published-minus-holdout AUC attenuation: {statistics.mean(atts):.6f}",
        f"Median holdout-minus-preholdout independent AUC delta: {statistics.median(deltas):.6f}",
        f"Holdout vs pre-holdout pooled-g sign concordance YES: {emp['YES']} / {len(avail)}",
        f"Holdout vs pre-holdout pooled-g sign concordance NO: {emp['NO']} / {len(avail)}",
        f"Holdout vs pre-holdout pooled-g sign concordance TIE: {emp['TIE']} / {len(avail)}",
    ]
summary += [
    "",
    "DEDUPLICATED GENE×ORIENTATION RESULTS",
    f"Available unique gene×orientation metrics: {len(unique_avail)} / 36",
    f"Directional AUC >0.5: {uniq_buckets['GT_0.5']} / {len(unique_avail)}",
    f"Directional AUC =0.5: {uniq_buckets['EQ_0.5']} / {len(unique_avail)}",
    f"Directional AUC <0.5: {uniq_buckets['LT_0.5']} / {len(unique_avail)}",
]
if unique_avail:
    uaucs=[float(r["DirectionalAUC"]) for r in unique_avail]
    ugs=[float(r["HedgesG_Aligned"]) for r in unique_avail]
    summary += [
        f"Median unique directional AUC: {statistics.median(uaucs):.6f}",
        f"Unique AUC range: {min(uaucs):.6f} to {max(uaucs):.6f}",
        f"Median unique aligned Hedges g: {statistics.median(ugs):.6f}",
    ]
summary += [
    "",
    "SPECIAL PRESPECIFIED REPORTING",
    tl_line,
    "",
    "UNAVAILABLE TARGETS",
]
if unavail:
    for r in unavail:
        summary.append(f"  PMID {r['PMID']} / {r['GeneSymbol']}: {r['Availability']}")
else:
    summary.append("  NONE")
summary += [
    "",
    "INTERPRETIVE GATE",
    "This is the standalone final temporal holdout. Stage4D remains the frozen development synthesis.",
    "Any later synthesis including GSE295540 must be a separately labeled post-holdout secondary analysis.",
]
(OUT/"STAGE5B_SUMMARY.txt").write_text("\n".join(summary)+"\n",encoding="utf-8")

# Freeze exact holdout outputs before any later updated synthesis.
manifest=[]
for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {
        "STAGE5B_FILE_MANIFEST.csv",
        "Stage5B_FinalTemporalHoldout_Bundle.zip",
        "STAGE5B_COMPLETE.flag"
    }:
        manifest.append({
            "RelativePath":str(p.relative_to(OUT)),
            "Bytes":p.stat().st_size,
            "SHA256":sha256(p)
        })
wcsv(OUT/"STAGE5B_FILE_MANIFEST.csv",manifest)

bundle=OUT/"Stage5B_FinalTemporalHoldout_Bundle.zip"
if bundle.exists(): bundle.unlink()
with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p!=bundle:
            z.write(p,p.relative_to(OUT))

COMPLETE.write_text(
    f"Completed {datetime.now(timezone.utc).isoformat()}. "
    "Standalone final temporal holdout frozen. Stage4D not updated.\n",
    encoding="utf-8"
)

print("\n==========================================================")
print("STAGE 5B COMPLETE")
print("==========================================================")
print("\n".join(summary))
print("\nResults:")
print(OUT/"GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv")
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
    if($LASTEXITCODE -ne 0){throw 'Stage 5B syntax check failed before unseal.'}
    & $PyExe -3 $Py
}else{
    & $PyExe -m py_compile $Py
    if($LASTEXITCODE -ne 0){throw 'Stage 5B syntax check failed before unseal.'}
    & $PyExe $Py
}
if($LASTEXITCODE -ne 0){throw "Stage 5B failed with exit code $LASTEXITCODE"}

Stop-Transcript
