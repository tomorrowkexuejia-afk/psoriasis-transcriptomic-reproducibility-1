from pathlib import Path
import csv, json, sys, math
from collections import defaultdict, Counter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from recompute import read_csv, write_csv, metrics, corrected_meta, classification

OUT = ROOT/'audit_impact'
OUT.mkdir(exist_ok=True)

# Inputs
pre_meta = read_csv(ROOT/'tables/S5_original_grid_corrected_meta_NOT_INDEPENDENCE_CLEARED.csv')
post_meta = read_csv(ROOT/'tables/source_adjudicated_subset_development_meta.csv')
post_inputs = read_csv(ROOT/'tables/source_adjudicated_subset_cohort_inputs.csv')
disposition = read_csv(ROOT/'tables/source_adjudicated_subset_claim_disposition.csv')
grid = read_csv(ROOT/'inputs/upload/planned_single_gene_tests_STAGE4_INPUT_FROZEN.csv')
mat544 = read_csv(ROOT/'inputs/data_review/matrices/GSE54456_primary_gene_matrix_STAGE3C.csv')

# helpers
key = lambda r: (str(r['PMID']), r['Gene'])
def fl(x): return float(x)
def cls(lo, hi): return classification(float(lo), float(hi))
def pct(x): return 100.0*x

# ---- 1. Provenance/source-use paired impact ----
pre_map = {key(r): r for r in pre_meta}
post_map = {key(r): r for r in post_meta}
source_rows=[]
for k, after in sorted(post_map.items()):
    before=pre_map[k]
    row={
        'PMID':k[0], 'Gene':k[1],
        'K_Before':int(before['K']), 'K_After':int(after['K']),
        'PooledG_Before':fl(before['HedgesG']), 'PooledG_After':fl(after['HedgesG']),
        'DeltaG':fl(after['HedgesG'])-fl(before['HedgesG']),
        'NormalCI_Before':cls(before['CI_Lower'], before['CI_Upper']),
        'NormalCI_After':cls(after['CI_Lower'], after['CI_Upper']),
        'mHK_Before':cls(before['mHK_Lower'], before['mHK_Upper']),
        'mHK_After':cls(after['mHK_Lower'], after['mHK_Upper']),
        'CIWidth_Before':fl(before['CI_Upper'])-fl(before['CI_Lower']),
        'CIWidth_After':fl(after['CI_Upper'])-fl(after['CI_Lower']),
    }
    row['ChangedByProvenanceAudit'] = abs(row['DeltaG']) > 1e-12 or row['K_Before'] != row['K_After']
    source_rows.append(row)
write_csv(OUT/'source_use_paired_33_claims.csv', source_rows)

not_eval=[]
for r in disposition:
    if str(r['IncludedInSynthesis']).lower() in ('false','0','no'):
        before=pre_map.get((str(r['PMID']),r['Gene']))
        not_eval.append({
            'PMID':str(r['PMID']), 'Gene':r['Gene'],
            'CandidateK':int(before['K']) if before else '',
            'RemainingEligibleCohorts':int(r['RemainingCohorts']),
            'Status':r['Status'],
            'Interpretation':'NOT_SUFFICIENTLY_EVALUABLE' if r['Status']!='PENDING_SUPPLEMENT' else 'PROVENANCE_UNRESOLVED'
        })
write_csv(OUT/'provenance_not_sufficiently_evaluable.csv', not_eval)

# ---- 2. Sample-overlap paired impact on fixed final cohort/claim set ----
orient={}
for r in grid:
    orient[(str(r['PMID']), r['GeneSymbol'])]=r['FrozenOrientation']

# full GSE54456 metrics, without the 42-sample exclusion, for exactly the final eligible GSE54456 contributions
full544={}
for r in post_inputs:
    if r['GSE']!='GSE54456': continue
    k=(str(r['PMID']),r['Gene'])
    sign=1 if orient[k]=='UP_IN_PSORIASIS' else -1
    cases=[sign*float(x[r['Gene']]) for x in mat544 if x['FinalPhenotype'].startswith('PSORIASIS_LESIONAL_SKIN')]
    ctrls=[sign*float(x[r['Gene']]) for x in mat544 if x['FinalPhenotype']=='HEALTHY_CONTROL_SKIN']
    full544[k]=metrics(cases,ctrls)

# group final corrected cohort inputs; create uncorrected copy by replacing only GSE54456 metrics
byclaim=defaultdict(list)
for r in post_inputs: byclaim[key(r)].append(r)
overlap_rows=[]
for k, corrected_rows in sorted(byclaim.items()):
    if k not in full544: continue
    uncorrected=[]
    for r in corrected_rows:
        if r['GSE']!='GSE54456':
            uncorrected.append({'GSE':r['GSE'],'HedgesG_Aligned':fl(r['HedgesG_Aligned']),'HedgesG_SE':fl(r['HedgesG_SE'])})
        else:
            m=full544[k]
            uncorrected.append({'GSE':'GSE54456','HedgesG_Aligned':m['HedgesG'],'HedgesG_SE':m['HedgesG_SE']})
    corrected=[{'GSE':r['GSE'],'HedgesG_Aligned':fl(r['HedgesG_Aligned']),'HedgesG_SE':fl(r['HedgesG_SE'])} for r in corrected_rows]
    mu=corrected_meta(uncorrected); mc=corrected_meta(corrected)
    def wshare(rows,m):
        w=np.array([1/(float(r['HedgesG_SE'])**2+float(m['Tau2'])) for r in rows])
        w=w/w.sum(); idx=[r['GSE'] for r in rows].index('GSE54456'); return float(w[idx])
    corr544=next(r for r in corrected_rows if r['GSE']=='GSE54456')
    unc544=full544[k]
    overlap_rows.append({
        'PMID':k[0], 'Gene':k[1], 'K':len(corrected_rows),
        'GSE54456_N_Before':unc544['N_Psoriasis']+unc544['N_Healthy'],
        'GSE54456_N_After':int(corr544['N_Psoriasis'])+int(corr544['N_Healthy']),
        'GSE54456_G_Before':unc544['HedgesG'],'GSE54456_G_After':fl(corr544['HedgesG']),
        'GSE54456_SE_Before':unc544['HedgesG_SE'],'GSE54456_SE_After':fl(corr544['HedgesG_SE']),
        'GSE54456_REWeightShare_Before':wshare(uncorrected,mu),
        'GSE54456_REWeightShare_After':wshare(corrected,mc),
        'PooledG_Before':mu['HedgesG'],'PooledG_After':mc['HedgesG'],
        'DeltaPooledG':mc['HedgesG']-mu['HedgesG'],
        'SE_Before':mu['SE'],'SE_After':mc['SE'],
        'CIWidth_Before':mu['CI_Upper']-mu['CI_Lower'],
        'CIWidth_After':mc['CI_Upper']-mc['CI_Lower'],
        'NormalCI_Before':classification(mu['CI_Lower'],mu['CI_Upper']),
        'NormalCI_After':classification(mc['CI_Lower'],mc['CI_Upper']),
        'mHK_Before':classification(mu['mHK_Lower'],mu['mHK_Upper']),
        'mHK_After':classification(mc['mHK_Lower'],mc['mHK_Upper']),
        'NormalPI_Before':classification(mu['NormalPI_Lower'],mu['NormalPI_Upper']),
        'NormalPI_After':classification(mc['NormalPI_Lower'],mc['NormalPI_Upper']),
        'tPI_Before':classification(mu['mHK_t_PI_Lower'],mu['mHK_t_PI_Upper']),
        'tPI_After':classification(mc['mHK_t_PI_Lower'],mc['mHK_t_PI_Upper']),
    })
write_csv(OUT/'overlap_correction_paired_30_claims.csv', overlap_rows)

# ---- 3. Direction-freezing illustrative counterfactual ----
direction_rows=[]
for r in post_inputs:
    auc=fl(r['AUC']); flipped=max(auc,1-auc)
    direction_rows.append({
        'PMID':str(r['PMID']),'Gene':r['Gene'],'GSE':r['GSE'],
        'FrozenDirectionalAUC':auc,'PostHocMaxAUC':flipped,
        'WouldBeFlipped':auc<0.5,
        'HedgesG_Aligned':fl(r['HedgesG_Aligned'])
    })
write_csv(OUT/'direction_freezing_counterfactual_176_contributions.csv', direction_rows)

# ---- 4. Inference-layer paired classifications ----
inference_rows=[]
for r in post_meta:
    inference_rows.append({
        'PMID':str(r['PMID']),'Gene':r['Gene'],
        'NormalCI':cls(r['CI_Lower'],r['CI_Upper']),
        'mHKCI':cls(r['mHK_Lower'],r['mHK_Upper']),
        'NormalPI':cls(r['NormalPI_Lower'],r['NormalPI_Upper']),
        'tBasedPI':cls(r['mHK_t_PI_Lower'],r['mHK_t_PI_Upper'])
    })
write_csv(OUT/'inference_layer_paired_33_claims.csv', inference_rows)

# summaries
src_changed=[r for r in source_rows if r['ChangedByProvenanceAudit']]
ov=np.array([abs(r['DeltaPooledG']) for r in overlap_rows])
ciabs=np.array([abs(r['CIWidth_After']-r['CIWidth_Before']) for r in overlap_rows])
w_before=np.array([r['GSE54456_REWeightShare_Before'] for r in overlap_rows]); w_after=np.array([r['GSE54456_REWeightShare_After'] for r in overlap_rows])
flipped=[r for r in direction_rows if r['WouldBeFlipped']]
inf_counts={}
for col in ['NormalCI','mHKCI','NormalPI','tBasedPI']:
    inf_counts[col]=dict(Counter(r[col] for r in inference_rows))

def n_changes(a,b): return sum(r[a]!=r[b] for r in inference_rows)
summary={
 'ProvenanceAudit':{
   'SourceUsedCohortTestsExcluded':17,
   'PairedFinalClaims':len(source_rows),
   'PairedClaimsWithChangedCohortSet':len(src_changed),
   'PairedClaimsUnchanged':len(source_rows)-len(src_changed),
   'ChangedClaims':[{'PMID':r['PMID'],'Gene':r['Gene'],'K_Before':r['K_Before'],'K_After':r['K_After'],'G_Before':r['PooledG_Before'],'G_After':r['PooledG_After'],'NormalCI_Before':r['NormalCI_Before'],'NormalCI_After':r['NormalCI_After'],'mHK_Before':r['mHK_Before'],'mHK_After':r['mHK_After']} for r in src_changed],
   'NotSufficientlyEvaluable':[r for r in not_eval]
 },
 'OverlapCorrection':{
   'AffectedClaims':len(overlap_rows),
   'GSE54456SampleN_Before':174,'GSE54456SampleN_After':132,'RemovedN':42,'SampleReductionPercent':100*(174-132)/174,
   'MedianAbsDeltaPooledG':float(np.median(ov)),'MaxAbsDeltaPooledG':float(np.max(ov)),
   'MedianAbsDeltaNormalCIWidth':float(np.median(ciabs)),
   'MedianGSE54456REWeightShare_Before':float(np.median(w_before)),
   'MedianGSE54456REWeightShare_After':float(np.median(w_after)),
   'NormalCIClassificationChanges':sum(r['NormalCI_Before']!=r['NormalCI_After'] for r in overlap_rows),
   'mHKClassificationChanges':sum(r['mHK_Before']!=r['mHK_After'] for r in overlap_rows),
   'NormalPIClassificationChanges':sum(r['NormalPI_Before']!=r['NormalPI_After'] for r in overlap_rows),
   'tPIClassificationChanges':sum(r['tPI_Before']!=r['tPI_After'] for r in overlap_rows)
 },
 'DirectionFreezing':{
   'Contributions':len(direction_rows),'BelowChanceContributions':len(flipped),
   'AffectedClaims':len(set((r['PMID'],r['Gene']) for r in flipped)),
   'AffectedGenes':sorted(set(r['Gene'] for r in flipped)),
   'FrozenAUCMedian':float(np.median([r['FrozenDirectionalAUC'] for r in direction_rows])),
   'PostHocMaxAUCMedian':float(np.median([r['PostHocMaxAUC'] for r in direction_rows])),
   'FrozenAUCRangeAffected':[float(min(r['FrozenDirectionalAUC'] for r in flipped)),float(max(r['FrozenDirectionalAUC'] for r in flipped))],
   'PostHocAUCRangeAffected':[float(min(r['PostHocMaxAUC'] for r in flipped)),float(max(r['PostHocMaxAUC'] for r in flipped))]
 },
 'InferenceLayer':{
   'Counts':inf_counts,
   'NormalCI_to_mHK_Changes':n_changes('NormalCI','mHKCI'),
   'NormalCI_to_NormalPI_Changes':n_changes('NormalCI','NormalPI'),
   'NormalPI_to_tPI_Changes':n_changes('NormalPI','tBasedPI'),
   'NormalCI_to_mHK_ChangedClaims':[{'PMID':r['PMID'],'Gene':r['Gene'],'Before':r['NormalCI'],'After':r['mHKCI']} for r in inference_rows if r['NormalCI']!=r['mHKCI']],
   'NormalCI_to_NormalPI_ChangedClaims':[{'PMID':r['PMID'],'Gene':r['Gene'],'Before':r['NormalCI'],'After':r['NormalPI']} for r in inference_rows if r['NormalCI']!=r['NormalPI']],
   'NormalPI_to_tPI_ChangedClaims':[{'PMID':r['PMID'],'Gene':r['Gene'],'Before':r['NormalPI'],'After':r['tBasedPI']} for r in inference_rows if r['NormalPI']!=r['tBasedPI']]
 }
}
(OUT/'audit_impact_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

# main-table rows, designed for manuscript insertion
changed=summary['ProvenanceAudit']['ChangedClaims']
main=[
 {'Audit component':'Provenance / source-use audit',
  'Paired comparison':'Original corrected candidate grid vs source-adjudicated evidence; final 33 claims compared on the same claim set',
  'Quantitative impact':'17 source-used cohort tests removed. Five PANoptosis claims fell from 4 candidate cohorts to 1 eligible cohort and became not sufficiently evaluable; INSIG1 was held because supplementary provenance remained unresolved. Among the 33 final claims, 31 were unchanged; GYS1 and TLN1 lost GSE121212. GYS1 pooled g 2.302 to 2.140; TLN1 -1.379 to -1.177. Normal-CI classifications did not change; TLN1 changed from negative to crossing zero under mHK.',
  'Interpretation':'The main effect was on independence and evaluability, not a broad loss of directional signal.'},
 {'Audit component':'Cross-accession sample-overlap correction',
  'Paired comparison':'Fixed final claim/cohort set; full GSE54456 (92 psoriasis, 82 healthy) vs removal of 42 samples shared with GSE13355',
  'Quantitative impact':f"30 claims affected. GSE54456 n fell from 174 to 132 (-24.1%). Median absolute change in pooled g was {summary['OverlapCorrection']['MedianAbsDeltaPooledG']:.3f} (maximum {summary['OverlapCorrection']['MaxAbsDeltaPooledG']:.3f}); median absolute change in normal-CI width was {summary['OverlapCorrection']['MedianAbsDeltaNormalCIWidth']:.3f}. Median random-effects weight share for GSE54456 changed from {100*summary['OverlapCorrection']['MedianGSE54456REWeightShare_Before']:.2f}% to {100*summary['OverlapCorrection']['MedianGSE54456REWeightShare_After']:.2f}%. No normal-CI, mHK, normal-PI, or t-based-PI classification changed.",
  'Interpretation':'Overlap correction changed the independence structure and some weights/effect estimates, while overall inferential classifications remained stable.'},
 {'Audit component':'Direction freezing',
  'Paired comparison':'Same 176 final cohort contributions; frozen directional AUC vs illustrative post hoc max(AUC, 1-AUC)',
  'Quantitative impact':f"{len(flipped)}/176 contributions ({100*len(flipped)/len(direction_rows):.1f}%), all from TLN1, had directional AUC <0.5. Retrospective reversal would convert AUCs {summary['DirectionFreezing']['FrozenAUCRangeAffected'][0]:.3f}-{summary['DirectionFreezing']['FrozenAUCRangeAffected'][1]:.3f} to {summary['DirectionFreezing']['PostHocAUCRangeAffected'][0]:.3f}-{summary['DirectionFreezing']['PostHocAUCRangeAffected'][1]:.3f}. The overall median AUC changed only from {summary['DirectionFreezing']['FrozenAUCMedian']:.3f} to {summary['DirectionFreezing']['PostHocMaxAUCMedian']:.3f}.",
  'Interpretation':'A global summary can look almost unchanged while a biologically important directional conflict is hidden.'},
 {'Audit component':'Inferential layer',
  'Paired comparison':'Same 33 final development claims under conventional CI, mHK CI, normal PI, and t-based PI',
  'Quantitative impact':'Normal CI: 32 positive, 1 negative. mHK CI: 30 positive, 3 crossing zero (3 classifications changed vs normal CI). Normal PI: 30 positive, 3 crossing zero (3 changed vs normal CI). t-based PI: 16 positive, 17 crossing zero (14 additional changes vs normal PI).',
  'Interpretation':'Direction was more stable than inferential certainty, especially when uncertainty for a new cohort was emphasized.'}
]
write_csv(OUT/'main_table4_audit_impact.csv',main)
print(json.dumps(summary,indent=2,ensure_ascii=False))
