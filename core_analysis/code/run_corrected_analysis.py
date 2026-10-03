"""Complete supplied-matrix arithmetic review. Original inputs are read-only.
Run with --project directory holding upload/ and data_review/matrices/.
"""
import argparse, csv, json, hashlib, sys, subprocess
from pathlib import Path
from collections import Counter,defaultdict
from recompute import read_csv,read_zip,write_csv,metrics,corrected_meta,classification

p=argparse.ArgumentParser()
p.add_argument('--project',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args(); root=a.project.resolve(); out=a.out.resolve(); out.mkdir(parents=True,exist_ok=True)
u=root/'upload'; mat=root/'data_review'/'matrices'
exfile=u/'GSE54456_STAGE4_EXCLUSION_SET_SHARED_42_FROZEN.csv'
gridfile=u/'planned_single_gene_tests_STAGE4_INPUT_FROZEN.csv'
mapfile=u/'GSE13355_GSE54456_DETERMINISTIC_SHARED_FROZEN.csv'
exc=read_csv(exfile); grid=read_csv(gridfile); mapping=read_csv(mapfile)
old=read_zip(u/'Stage4A_SingleGene_PerCohort_Bundle.zip','single_gene_per_cohort_metrics_STAGE4A.csv')
preold=read_zip(u/'Stage4F_ManuscriptReady_Bundle.zip','Table1_Primary_Marker_Results_STAGE4F.csv')
mats={x.name.split('_')[0]:read_csv(x) for x in mat.glob('*_primary_gene_matrix_STAGE3C.csv')}
ex={r['GSM'] for r in exc}; assert len(exc)==len(ex)==42
by={g:{r['GSM']:r for r in rs} for g,rs in mats.items()}
assert len(by['GSE13355'])==len(mats['GSE13355']) and len(by['GSE54456'])==len(mats['GSE54456'])
for r in exc:
 s=by['GSE54456'][r['GSM']]
 assert r['FinalPhenotype']==s['FinalPhenotype'] and r['OfficialOverlapToken']==s['Title']
for r in mapping:
 for g in ['GSE13355','GSE54456']:
  s=by[g][r[g+'_GSM']]
  assert r[g+'_Phenotype']==s['FinalPhenotype'] and r[g+'_Title']==s['Title']
 assert r['GSE13355_Phenotype']==r['GSE54456_Phenotype']
for g in ['GSE13355','GSE54456']:
 assert len({r[g+'_GSM'] for r in mapping})==len(mapping)
matched={r['GSE54456_GSM'] for r in mapping}; assert matched<=ex
key=lambda r:tuple(r[c] for c in ['PMID','GeneSymbol','FrozenOrientation','GSE','SampleExclusionSetID'])
assert Counter(key(r) for r in grid)==Counter(key(r) for r in old)
lookup={key(r):r for r in old}
assert all(int(r['RequiredExcludedN'])==int(lookup[key(r)]['ExcludedN']) for r in grid)
write_csv(out/'overlap_sample_check.csv',[dict(r,DirectCrossCohortMatchSupplied=r['GSM'] in matched) for r in exc])
subprocess.run([sys.executable,str(Path(__file__).with_name('recompute.py')),'--matrices',str(mat),'--stage4a',str(u/'Stage4A_SingleGene_PerCohort_Bundle.zip'),'--stage4f',str(u/'Stage4F_ManuscriptReady_Bundle.zip'),'--stage5b',str(u/'Stage5B_FinalTemporalHoldout_Bundle.zip'),'--holdout',str(u/'GSE295540_raw_counts_All_samples.csv.gz'),'--exclusions',str(exfile),'--out',str(out)],check=True,stdout=subprocess.DEVNULL)
hold={(r['PMID'],r['Gene']):r for r in read_csv(out/'holdout_recomputed.csv')}
oldgroups=defaultdict(list)
for r in old: oldgroups[(r['PMID'],r['GeneSymbol'])].append(r)
inputs=[]; before=[]; after=[]; loco=[]; platform=[]
for (pmid,gene),rs in oldgroups.items():
 both={'GSE13355','GSE54456'}<={r['GSE'] for r in rs}
 sr=[]
 for r in rs:
  remove=r['GSE']=='GSE54456' and (both or int(r['ExcludedN'])>0)
  data=[x for x in mats[r['GSE']] if not remove or x['GSM'] not in ex]
  sign=1 if r['FrozenOrientation']=='UP_IN_PSORIASIS' else -1
  m=metrics([sign*float(x[gene]) for x in data if x['FinalPhenotype'].startswith('PSORIASIS_LESIONAL_SKIN')], [sign*float(x[gene]) for x in data if x['FinalPhenotype']=='HEALTHY_CONTROL_SKIN'])
  row=dict(PMID=pmid,Gene=gene,GSE=r['GSE'],FrozenOrientation=r['FrozenOrientation'],ExcludedN=42 if remove else 0,**m,HedgesG_Aligned=m['HedgesG'])
  sr.append(row); inputs.append(row)
 def labeled(rows):
  m=corrected_meta(rows)
  return dict(PMID=pmid,Gene=gene,**m,PrimaryClassification=classification(m['CI_Lower'],m['CI_Upper']),mHK_Classification=classification(m['mHK_Lower'],m['mHK_Upper']),NormalPIClassification=classification(m['NormalPI_Lower'],m['NormalPI_Upper']),PostHoc_t_PIClassification=classification(m['mHK_t_PI_Lower'],m['mHK_t_PI_Upper']))
 before.append(labeled(sr))
 hr=hold[(pmid,gene)]
 after.append(labeled(sr+[dict(HedgesG_Aligned=hr['HedgesG'],HedgesG_SE=hr['HedgesG_SE'])]))
 for removed in sr:
  q=labeled([r for r in sr if r['GSE']!=removed['GSE']]); loco.append(dict(RemovedGSE=removed['GSE'],**q))
 buckets=defaultdict(list)
 for r in sr: buckets['MICROARRAY' if r['GSE'] in ['GSE13355','GSE14905','GSE78097','GSE201827'] else 'RNA_SEQ'].append(r)
 for name,rows in buckets.items():
  if len(rows)>=2: platform.append(dict(Platform=name,**labeled(rows)))
write_csv(out/'synthesis_per_cohort_inputs_corrected.csv',inputs)
write_csv(out/'development_meta_all39_corrected.csv',before)
write_csv(out/'postholdout_meta_all39_corrected.csv',after)
write_csv(out/'development_leave_one_cohort_out_corrected.csv',loco)
write_csv(out/'development_platform_strata_corrected.csv',platform)
oldmap={(r['PMID'],r['GeneSymbol']):r for r in preold}
comp=[]
for r in before:
 o=oldmap[(r['PMID'],r['Gene'])]
 oc=classification(float(o['PooledHedgesG_CI95_Lower']),float(o['PooledHedgesG_CI95_Upper']))
 comp.append(dict(PMID=r['PMID'],Gene=r['Gene'],OldG=o['PooledHedgesG_REML'],CorrectedG=r['HedgesG'],OldTau2=o['Tau2_REML'],CorrectedTau2=r['Tau2'],PrimaryClassificationChanged=oc!=r['PrimaryClassification']))
write_csv(out/'development_meta_changes.csv',comp)
# Retain a conventional Stage5B-shaped corrected table for downstream use.
ho=read_zip(u/'Stage5B_FinalTemporalHoldout_Bundle.zip','GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv')
fields={'DirectionalAUC':'AUC','AUC_SE_DeLong':'AUC_SE','AUC_CI95_Lower':'AUC_Lower','AUC_CI95_Upper':'AUC_Upper','HedgesG_Aligned':'HedgesG','HedgesG_SE':'HedgesG_SE','HedgesG_CI95_Lower':'HedgesG_Lower','HedgesG_CI95_Upper':'HedgesG_Upper','OrientedPsoriasisMean':'OrientedPsoriasisMean','OrientedHealthyMean':'OrientedHealthyMean','OrientedMeanDifference':'OrientedMeanDifference'}
beforemap={(r['PMID'],r['Gene']):r for r in before}
for r in ho:
 h=hold[(r['PMID'],r['GeneSymbol'])]
 for dest,src in fields.items(): r[dest]=h[src]
 pre=beforemap[(r['PMID'],r['GeneSymbol'])]['HedgesG']; r['PreHoldoutPooledHedgesG']=pre
 r['DirectionConcordantWithPublished']='YES' if float(h['HedgesG'])>0 else 'NO'
 r['HoldoutVsPreHoldoutEmpiricalSignConcordance']='YES' if float(h['HedgesG'])*pre>0 else 'NO'
write_csv(out/'GSE295540_holdout_results_CORRECTED_REVIEW.csv',ho)
counts=lambda rs:{c:dict(Counter(r[c] for r in rs)) for c in ['PrimaryClassification','mHK_Classification','NormalPIClassification','PostHoc_t_PIClassification']}
loco_bad=[r for r in loco if r['PrimaryClassification']!=beforemap[(r['PMID'],r['Gene'])]['PrimaryClassification']]
pby=defaultdict(list)
for r in platform:pby[(r['PMID'],r['Gene'])].append(r)
paired=[rs for rs in pby.values() if len(rs)==2]
s=dict(Development=counts(before),PostHoldout=counts(after),DevelopmentPrimaryChanges=sum(r['PrimaryClassificationChanged'] for r in comp),LOCOTests=len(loco),LOCOPrimaryClassificationChanges=len(loco_bad),PlatformPairedGroups=len(paired),PlatformSignDiscordant=sum(rs[0]['HedgesG']*rs[1]['HedgesG']<0 for rs in paired),ExcludedSamples=len(ex),DirectMappedSamples=len(mapping),ExtraExcludedSamples=[r['GSM'] for r in exc if r['GSM'] not in matched],AllInputsRecomputedFromMatrices=True,Scope='Provided processed development matrices and raw-count holdout; this numeric module retains the original grid and does not clear independence. External preprocessing and provenance checks are separately recorded under evidence/; use amended tables for interpretation.')
(out/'complete_summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
# These former files have an explicitly different scope; remove stale partial-review outputs from this new review directory.
for name in ['posthoc_prediction_interval_sensitivity.csv','development_pending_exclusions.csv','reml_pending_groups.csv']:
 (out/name).unlink(missing_ok=True)
checks=read_csv(out/'input_checksums.csv')
for f in [gridfile,mapfile]:checks.append(dict(Filename=f.name,Bytes=f.stat().st_size,SHA256=hashlib.sha256(f.read_bytes()).hexdigest()))
write_csv(out/'input_checksums.csv',checks)
print(json.dumps(s,ensure_ascii=False,indent=2))
