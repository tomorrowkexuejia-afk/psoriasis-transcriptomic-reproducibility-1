"""Build evidence-amended analyses, manuscript tables and publication figures.
The amendment is post hoc. Historical numeric reproduction is retained separately.
"""
from pathlib import Path
import argparse,csv,json,statistics,math
from collections import defaultdict,Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from recompute import read_csv,write_csv,metrics,corrected_meta
p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);a=p.parse_args();root=a.release
rs=root/'results';out=root/'tables';fig=root/'figures';out.mkdir(exist_ok=True);fig.mkdir(exist_ok=True)
source=read_csv(rs/'development_recomputed.csv');grid=read_csv(root/'inputs/upload/planned_single_gene_tests_STAGE4_INPUT_FROZEN.csv')
old=defaultdict(list)
for r in grid:old[(r['PMID'],r['GeneSymbol'])].append(r)
mat={f.name.split('_')[0]:read_csv(f) for f in (root/'inputs/data_review/matrices').glob('*.csv')}
ex={r['GSM'] for r in read_csv(root/'inputs/upload/GSE54456_STAGE4_EXCLUSION_SET_SHARED_42_FROZEN.csv')}
hold={(r['PMID'],r['Gene']):r for r in read_csv(rs/'holdout_recomputed.csv')}
used={'39480805':{'GSE54456','GSE66511','GSE121212'},'41224720':{'GSE121212'}}
uncertain={'39735895'}
adjudication=[]
for r in grid:
 pm=r['PMID'];g=r['GSE'];bad=g in used.get(pm,set());pending=pm in uncertain and not bad
 status='EXCLUDE_SOURCE_USED' if bad else 'PENDING_SUPPLEMENT_ACCESSIONS' if pending else 'NO_ADDITIONAL_CONFLICT_IDENTIFIED'
 note='Source Table 1 and Methods include this validation cohort' if pm=='39480805' and bad else 'Source Figure 5 B/C explicitly labels GSE121212' if bad else 'Source reports additional validation cohorts; complete supplementary accession audit outstanding' if pending else 'Original exclusion rules retained; full-text and available figure accession review found no additional selected-cohort conflict'
 adjudication.append(dict(PMID=pm,Gene=r['GeneSymbol'],GSE=g,AmendmentStatus=status,Reason=note))
write_csv(out/'S1_eligibility_amendment.csv',adjudication)
phases={};cohort_tables={};claim_tables={};loco=[];platform=[];small=[]
for mode in ['confirmed_use_corrections','source_adjudicated_subset']:
 full=[];summ=[];post=[];claims=[]
 for (pm,gene),orig in old.items():
  selected=[r for r in orig if r['GSE'] not in used.get(pm,set()) and not(mode=='source_adjudicated_subset' and pm in uncertain)]
  eligible=len(selected)>=3
  claims.append(dict(PMID=pm,Gene=gene,RemainingCohorts=len(selected),IncludedInSynthesis=eligible,Status='INCLUDED' if eligible else 'PENDING_SUPPLEMENT' if mode=='source_adjudicated_subset' and pm in uncertain else 'FEWER_THAN_3_COHORTS'))
  if not eligible:continue
  both={'GSE13355','GSE54456'}<={r['GSE'] for r in selected};inputrows=[]
  for r in selected:
   sign=1 if r['FrozenOrientation']=='UP_IN_PSORIASIS' else -1
   remove=r['GSE']=='GSE54456' and (both or int(r['RequiredExcludedN']))
   data=[x for x in mat[r['GSE']] if not remove or x['GSM'] not in ex]
   m=metrics([sign*float(x[gene]) for x in data if x['FinalPhenotype'].startswith('PSORIASIS_LESIONAL_SKIN')],[sign*float(x[gene]) for x in data if x['FinalPhenotype']=='HEALTHY_CONTROL_SKIN'])
   q=dict(PMID=pm,Gene=gene,GSE=r['GSE'],RemovedSamples=42 if remove else 0,**m,HedgesG_Aligned=m['HedgesG']);inputrows.append(q);full.append(q)
  def meta(rows):
   m=corrected_meta(rows)
   return dict(PMID=pm,Gene=gene,**m)
  summ.append(meta(inputrows));h=hold[(pm,gene)];post.append(meta(inputrows+[dict(HedgesG_Aligned=h['HedgesG'],HedgesG_SE=h['HedgesG_SE'])]))
  if mode=='source_adjudicated_subset':
   for x in inputrows:loco.append(dict(RemovedGSE=x['GSE'],**meta([r for r in inputrows if r['GSE']!=x['GSE']])))
   for plat,gses in [('Microarray',{'GSE13355','GSE14905','GSE78097','GSE201827'}),('RNA-seq',{'GSE54456','GSE66511','GSE121212'})]:
    group=[r for r in inputrows if r['GSE'] in gses]
    if len(group)>=2:platform.append(dict(Platform=plat,**meta(group)))
   group=[r for r in inputrows if r['N_Healthy']>=10]
   if len(group)>=2:small.append(meta(group))
 def counts(items):
  return {n:dict(Counter('POSITIVE' if r[lo]>0 else 'NEGATIVE' if r[hi]<0 else 'CROSSES_ZERO' for r in items)) for n,lo,hi in [('NormalCI','CI_Lower','CI_Upper'),('mHK','mHK_Lower','mHK_Upper'),('NormalPI','NormalPI_Lower','NormalPI_Upper'),('PostHoc_t_PI','mHK_t_PI_Lower','mHK_t_PI_Upper')]}
 hrows=[hold[(r['PMID'],r['Gene'])] for r in summ]
 phases[mode]=dict(Claims=len(summ),Genes=len({r['Gene'] for r in summ}),Papers=len({r['PMID'] for r in summ}),Rows=len(full),Development=counts(summ),PostHoldout=counts(post),MedianG=statistics.median(r['HedgesG'] for r in summ),MedianI2=statistics.median(r['I2'] for r in summ),MedianHoldoutAUC=statistics.median(float(r['AUC']) for r in hrows),MedianHoldoutG=statistics.median(float(r['HedgesG']) for r in hrows),PublishedDirectionHoldout=sum(float(r['HedgesG'])>0 for r in hrows))
 write_csv(out/f'{mode}_cohort_inputs.csv',full);write_csv(out/f'{mode}_development_meta.csv',summ);write_csv(out/f'{mode}_postholdout_meta.csv',post);write_csv(out/f'{mode}_claim_disposition.csv',claims)
 cohort_tables[mode]=full;claim_tables[mode]=summ
write_csv(out/'S7_adjudicated_leave_one_cohort_out.csv',loco);write_csv(out/'S8_adjudicated_platform.csv',platform);write_csv(out/'S9_adjudicated_small_control_exclusion.csv',small)
write_csv(out/'S2_original_202_tests_arithmetic_reproduction.csv',source)
write_csv(out/'S3_corrected_holdout_39_claims.csv',list(hold.values()))
unique={r['Gene']:r for r in hold.values()};write_csv(out/'S4_corrected_holdout_36_genes.csv',list(unique.values()))
# Ref-compatible full corrected table for the original 39-claim cohort grid is retained as a numeric audit only.
write_csv(out/'S5_original_grid_corrected_meta_NOT_INDEPENDENCE_CLEARED.csv',read_csv(rs/'development_meta_all39_corrected.csv'))
(root/'amended_summary.json').write_text(json.dumps(phases,indent=2))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(f,n):
 for ext in ['png','pdf','svg']:f.savefig(fig/(n+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(f)
# Figure 1: factual audit flow; no fabricated exclusion reasons for original screening.
f,ax=plt.subplots(figsize=(8.5,6));ax.axis('off')
boxes=[(.05,.75,.9,.17,'Candidate inventory: 15 papers / 39 claims / 36 genes\n202 candidate paper-gene-cohort tests'),(.05,.49,.43,.17,'Confirmed source reuse\n17 tests excluded\n5 claims have <3 cohorts'),(.53,.49,.42,.17,'Additional supplement audit\nINSIG1 supplementary accessions\n1 claim held pending'),(.05,.17,.9,.20,'Source-adjudicated subset: 13 papers / 33 claims / 32 genes\n176 cohort contributions; overlap removal applied within each synthesis\nSeparate RNA-seq holdout: 7 psoriasis / 6 healthy')]
for x,y,w,h,t in boxes:
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008',fc='#f3f5f7',ec='#415366'));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=10)
for xx in [.27,.74]:
 ax.annotate('',xy=(xx,.67),xytext=(.5,.75),arrowprops={'arrowstyle':'->','color':'#415366'});ax.annotate('',xy=(.5,.37),xytext=(xx,.49),arrowprops={'arrowstyle':'->','color':'#415366'})
ax.text(.5,.03,'Cohort eligibility follows source use; original score orientations are preserved.',ha='center',fontsize=9)
save(f,'Figure1_Provenance_Audit_Flow')
# Figure 2: corrected independent-subset directional AUC distributions.
cs=cohort_tables['source_adjudicated_subset'];gs=sorted({r['GSE'] for r in cs});f,ax=plt.subplots(figsize=(8.8,4.5));rng=np.random.default_rng(20260907)
vals=[[r['AUC'] for r in cs if r['GSE']==g] for g in gs]
ax.boxplot(vals,positions=np.arange(len(gs)),widths=.5,showfliers=False,patch_artist=True,boxprops={'facecolor':'#d8e4ec'},medianprops={'color':'black'})
for i,v in enumerate(vals):ax.scatter(i+rng.uniform(-.16,.16,len(v)),v,s=13,alpha=.65,color='#285a78')
ax.axhline(.5,ls='--',lw=1,color='#777');ax.set_xticks(range(len(gs)),[f'{g}\n(n={len(v)})' for g,v in zip(gs,vals)],fontsize=8);ax.set_ylabel('Directional AUC');ax.set_ylim(-.03,1.04);ax.set_title('Source-adjudicated subset: 176 cohort contributions',loc='left',fontweight='bold');save(f,'Figure2_Adjudicated_AUC')
# Figure 3: distinguish CIs and PIs; use all phase counts automatically.
f,axs=plt.subplots(1,2,figsize=(9,4.5));ph=phases['source_adjudicated_subset']
for ax,names,title in [(axs[0],[('Development','NormalCI'),('Development','mHK'),('PostHoldout','NormalCI'),('PostHoldout','mHK')],'Confidence intervals'),(axs[1],[('Development','NormalPI'),('Development','PostHoc_t_PI'),('PostHoldout','NormalPI'),('PostHoldout','PostHoc_t_PI')],'Prediction intervals')]:
 for i,(stage,metric) in enumerate(names):
  c=ph[stage][metric];pos=c.get('POSITIVE',0);cross=c.get('CROSSES_ZERO',0);neg=c.get('NEGATIVE',0)
  ax.barh(i,pos,color='#356e8d');ax.barh(i,cross,left=pos,color='#b7bcc1');ax.barh(i,neg,left=pos+cross,color='#ad5348')
  if pos:ax.text(pos/2,i,str(pos),va='center',ha='center',color='white')
  if cross:ax.text(pos+cross/2,i,str(cross),va='center',ha='center',fontsize=9)
  if neg:ax.text(pos+cross+neg/2,i,str(neg),va='center',ha='center',color='white',fontsize=8)
 ax.set_yticks(range(4),[('Dev.' if st=='Development' else '+ holdout')+' '+{'NormalCI':'normal','mHK':'mHK','NormalPI':'normal','PostHoc_t_PI':'t sensitivity'}[m] for st,m in names],fontsize=9);ax.invert_yaxis();ax.set_xlim(0,33);ax.set_xlabel('Paper–gene claims');ax.set_title(title,fontweight='bold')
f.text(.5,-.015,'Blue: interval entirely positive. Grey: crosses zero. Red: entirely negative. t-based PI: sensitivity analysis.',ha='center',fontsize=9);f.tight_layout(w_pad=2);save(f,'Figure3_Adjudicated_Interval_Sensitivity')
# Figure 4: TLN1 remains an explicitly provisional case study, not an independence-cleared primary finding.
f,axs=plt.subplots(1,2,figsize=(9,5),sharey=True)
for ax,gene in zip(axs,['TLN1','GYS1']):
 rows=[r for r in cs if r['Gene']==gene and r['PMID']=='41224720']+[r for r in source if r['Gene']==gene and r['PMID']=='41224720' and r['GSE']=='GSE121212'];rows=sorted(rows,key=lambda r:r['GSE'])
 for i,r in enumerate(rows):
  yy=float(r['HedgesG']);lo=float(r['HedgesG_Lower']);hi=float(r['HedgesG_Upper']);color='#ba4b42' if r['GSE']=='GSE121212' else '#777'
  ax.errorbar(yy,i,xerr=[[yy-lo],[hi-yy]],fmt='o',color=color,capsize=2)
 h=hold[('41224720',gene)];yy=float(h['HedgesG']);lo=float(h['HedgesG_Lower']);hi=float(h['HedgesG_Upper']);ax.errorbar(yy,7,xerr=[[yy-lo],[hi-yy]],fmt='s',color='#24688b',capsize=3)
 ax.axvline(0,color='black',lw=.8);ax.set_title(gene,fontweight='bold');ax.set_xlabel('Aligned Hedges g (95% CI)')
axs[0].set_yticks(range(8),[r['GSE'] for r in rows]+['GSE295540 holdout'],fontsize=9);axs[0].invert_yaxis();f.text(.5,.012,'Red: confirmed source reuse. Grey: eligible development cohort. Blue: separate holdout.\nTLN1 direction follows the source main text; its supplementary figure shows an opposing pattern.',ha='center',fontsize=9);f.tight_layout(rect=(0,.10,1,1));save(f,'Figure4_TLN1_GYS1_Provenance_Qualified')
# Supplementary forest, 33 claims.
r=sorted(claim_tables['source_adjudicated_subset'],key=lambda x:x['HedgesG']);f,ax=plt.subplots(figsize=(8,11))
for i,x in enumerate(r):ax.errorbar(x['HedgesG'],i,xerr=[[x['HedgesG']-x['CI_Lower']],[x['CI_Upper']-x['HedgesG']]],fmt='o',color='#356e8d',ms=4,capsize=2)
ax.set_yticks(range(len(r)),[f"{x['Gene']} | PMID {x['PMID']}" for x in r],fontsize=8);ax.axvline(0,color='black',lw=.8);ax.set_xlabel('Pooled aligned Hedges g (normal 95% CI)');ax.set_title('Source-adjudicated subset: 33 claims',loc='left');save(f,'FigureS1_Adjudicated_Forest')
print(json.dumps(phases,indent=2))
