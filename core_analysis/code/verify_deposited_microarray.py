import csv,gzip
from collections import defaultdict
from pathlib import Path
import argparse
p=argparse.ArgumentParser();p.add_argument("--deposits",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];DEPOSITS=a.deposits;OUTPUT=a.output;OUTPUT.mkdir(parents=True,exist_ok=True)
import numpy as np
nom=list(csv.DictReader(open(str(ROOT/'inputs/upload')+'/GSE295540_target_nomenclature_FINAL_FROZEN.csv',encoding='utf-8-sig')));targets={r['NCBICurrentSymbol'] for r in nom}
probes=defaultdict(set)
with gzip.open(str(DEPOSITS)+'/GPL570.annot.gz','rt') as f:
 for line in f:
  if line.startswith('!platform_table_begin'):break
 rd=csv.DictReader(f,delimiter='\t')
 for r in rd:
  if r['ID'].startswith('!platform_table_end'):break
  for symbol in r['Gene symbol'].split('///'):
   if symbol.strip() in targets:probes[r['ID']].add(symbol.strip())
results=[]
for g in ['GSE13355','GSE14905','GSE78097','GSE201827']:
 processed=list(csv.DictReader(open(str(ROOT/'inputs/data_review/matrices')+'/'+g+'_primary_gene_matrix_STAGE3C.csv',encoding='utf-8-sig')))
 values=defaultdict(list)
 with gzip.open(str(DEPOSITS)+'/'+g+'_series_matrix.txt.gz','rt') as f:
  for line in f:
   if line.startswith('!series_matrix_table_begin'):break
  rd=csv.reader(f,delimiter='\t');head=next(rd);indices=[head.index(r['GSM']) for r in processed]
  seen=set();complete=True
  try:
   for r in rd:
    if not r or r[0].startswith('!series_matrix_table_end'):break
    if r[0] in probes:
     x=[float(r[i]) for i in indices];seen.add(r[0])
     for gene in probes[r[0]]:values[gene].append(x)
  except EOFError:
   complete=False
  print(g, 'complete_table',complete,'target_probes_seen',len(seen),'expected',len(probes),flush=True)
 for gene in sorted(targets):
  if not values[gene]:results.append(dict(GSE=g,Gene=gene,Probes=0,MaxAbsoluteDifference='',Status='MAPPING_NOT_RECONSTRUCTED' if complete else 'INCOMPLETE_SOURCE_NOT_VERIFIED'));continue
  x=np.median(np.array(values[gene]),axis=0);d=max(abs(float(r[gene])-float(x[i])) for i,r in enumerate(processed))
  results.append(dict(GSE=g,Gene=gene,Probes=len(values[gene]),MaxAbsoluteDifference=d if complete else '',Status=('MATCH' if d<1e-6 else 'MISMATCH') if complete else 'INCOMPLETE_SOURCE_NOT_VERIFIED'))
 print(g,[(r['Gene'],r['MaxAbsoluteDifference']) for r in results if r['GSE']==g and r['Status']!='MATCH'],flush=True)
with open(str(OUTPUT/'deposited_microarray_preprocessing_check.csv'),'w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
