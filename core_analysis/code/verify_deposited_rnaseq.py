import csv,gzip,re,json,hashlib
from pathlib import Path
import argparse
p=argparse.ArgumentParser();p.add_argument("--deposits",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];DEPOSITS=a.deposits;OUTPUT=a.output;OUTPUT.mkdir(parents=True,exist_ok=True)
from collections import defaultdict
import numpy as np
out=OUTPUT; results=[]
nom=list(csv.DictReader(open(str(ROOT/'inputs/upload')+'/GSE295540_target_nomenclature_FINAL_FROZEN.csv',encoding='utf-8-sig')))
# Only the documented historical rename needed among the 36 target symbols in these deposits.
renames={'ARHGEF28':'RGNEF'}
for g,n,mode in [('GSE54456','GSE54456_RPKM_samples.txt.gz','RPKM'),('GSE66511','GSE66511_Psoriasis_counts.txt.gz','COUNTS'),('GSE121212','GSE121212_readcount.txt.gz','COUNTS')]:
 with gzip.open(str(DEPOSITS)+'/'+n,'rt') as f:
  rd=csv.reader(f,delimiter='\t');header=next(rd);raw=list(rd)
 val=np.array([r[1:] for r in raw],float);assert np.isfinite(val).all() and (val>=0).all()
 ids=defaultdict(list)
 for i,r in enumerate(raw):ids[r[0]].append(i)
 libs=val.sum(axis=0);processed=list(csv.DictReader(open(str(ROOT/'inputs/data_review/matrices')+'/'+g+'_primary_gene_matrix_STAGE3C.csv',encoding='utf-8-sig')))
 for r in nom:
  gene=r['NCBICurrentSymbol'];symbol=gene if gene in ids else renames.get(gene,gene)
  if symbol not in ids:results.append(dict(GSE=g,Gene=gene,RawSymbol=symbol,ComparedSamples=0,MaxAbsoluteDifference='',Status='MAPPING_NOT_RECONSTRUCTED'));continue
  z=val[ids[symbol]]
  x=np.median(np.log2(z+1),axis=0) if mode=='RPKM' else np.log2(z.sum(axis=0)/libs*1e6+1)
  differences=[abs(float(row[gene])-float(x[header.index(row['ExpressionColumn'])-1])) for row in processed]
  delta=max(differences);results.append(dict(GSE=g,Gene=gene,RawSymbol=symbol,ComparedSamples=len(differences),MaxAbsoluteDifference=delta,Status='MATCH' if delta<1e-6 else 'MISMATCH'))
 print(g,'maxdiff',max(float(r['MaxAbsoluteDifference']) for r in results if r['GSE']==g and r['MaxAbsoluteDifference']!=''),flush=True)
with open(out/'deposited_rnaseq_preprocessing_check.csv','w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
# Original overlap attachment directly supports all 42 exclusions, even when no one-to-one cross-map was supplied.
tokens=gzip.open(str(DEPOSITS)+'/GSE54456_MAoverlappedsamples.txt.gz','rt').read().split()
e=list(csv.DictReader(open(str(ROOT/'inputs/upload')+'/GSE54456_STAGE4_EXCLUSION_SET_SHARED_42_FROZEN.csv',encoding='utf-8-sig')))
assert len(tokens)==42 and set(tokens)=={r['OfficialOverlapToken'] for r in e}
(out/'official_overlap_check.json').write_text(json.dumps({'OfficialTokens':42,'FrozenExclusionTokens':42,'ExactSetMatch':True}))
