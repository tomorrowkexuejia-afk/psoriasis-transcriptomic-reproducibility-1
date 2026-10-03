"""Audit supplied frozen mappings; reproduce original alias-based holdout.
Usage: python audit_mapping.py --project PROJECT --out NEW_OUTPUT
Requires numpy and scipy, and recompute.py beside this file.
"""
import argparse,gzip,csv,json,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from recompute import read_csv,read_zip,write_csv,metrics
p=argparse.ArgumentParser();p.add_argument('--project',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
u=a.project/'upload';a.out.mkdir(parents=True,exist_ok=True)
nf=u/'GSE295540_target_nomenclature_FINAL_FROZEN.csv';sf=u/'GSE295540_sample_phenotype_map_FINAL_FROZEN.csv';raw=u/'GSE295540_raw_counts_All_samples.csv.gz';bundle=u/'Stage5B_FinalTemporalHoldout_Bundle.zip'
nom=read_csv(nf); sm=read_csv(sf); original=read_zip(bundle,'GSE295540_paper_gene_holdout_results_FINAL_STAGE5B.csv');oldmap=read_zip(bundle,'GSE295540_target_mapping_audit_STAGE5B.csv')
assert len(nom)==len({r['NCBIGeneID'] for r in nom})==len({r['NCBICurrentSymbol'] for r in nom})==36
assert len(sm)==len({r['GSM'] for r in sm})==len({r['FrozenLibraryAlias'] for r in sm})==13
assert Counter(r['FinalPhenotype'] for r in sm)=={'PSORIASIS_LESIONAL_SKIN':7,'HEALTHY_CONTROL_SKIN':6}
with gzip.open(raw,'rt',encoding='utf-8-sig') as f:
 rr=list(csv.reader(f));header=rr[0]; names=[r[0] for r in rr[1:]];v=np.array([r[1:] for r in rr[1:]],float)
assert header==['Gene']+[f'NL{i}' for i in range(1,7)]+[f'PSO{i}' for i in range(1,8)]
assert len(names)==len(set(names))==len(set(s.upper() for s in names))
assert np.isfinite(v).all() and (v>=0).all()
idx={s.upper():i for i,s in enumerate(names)}; libs=v.sum(axis=0); aliases=defaultdict(set)
for r in nom:
 for s in r['FrozenSymbolAliases'].split('|'):aliases[s.strip().upper()].add(r['NCBIGeneID'])
unique={s:next(iter(ids)) for s,ids in aliases.items() if len(ids)==1}
nomby={r['NCBIGeneID']:r for r in nom}; om={r['NCBIGeneID']:r for r in oldmap}; diag=[]; ori_metrics={}
for r in nom:
 gid=r['NCBIGeneID'];current=r['NCBICurrentSymbol'];assert current.upper() in idx
 hits=[i for i,n in enumerate(names) if unique.get(n.upper())==gid]
 assert len(hits)==int(om[gid]['RawRowsMapped'])
 x=np.log2(v[hits].sum(axis=0)/libs*1e6+1)
 direct=np.log2(v[idx[current.upper()]]/libs*1e6+1)
 ori_metrics[gid]=(x,direct)
 diag.append(dict(Gene=current,NCBIGeneID=gid,ExactCurrentSymbolRows=1,OriginalAliasMatchedRows=len(hits),OriginalMatchedSymbols='|'.join(names[i] for i in hits),ExtraSymbols='|'.join(names[i] for i in hits if names[i].upper()!=current.upper()),CorrectionNeeded=len(hits)!=1))
checks=[]
for r in original:
 gid=r['NCBIGeneID'];assert r['GeneSymbol']==nomby[gid]['NCBICurrentSymbol']
 assert r['FrozenOrientation'] in ['UP_IN_PSORIASIS','DOWN_IN_PSORIASIS']
 sign=1 if r['FrozenOrientation']=='UP_IN_PSORIASIS' else -1
 x,d=ori_metrics[gid]
 old=metrics(sign*x[6:],sign*x[:6]);direct=metrics(sign*d[6:],sign*d[:6])
 err=max(abs(old['AUC']-float(r['DirectionalAUC'])),abs(old['AUC_SE']-float(r['AUC_SE_DeLong'])),abs(old['HedgesG']-float(r['HedgesG_Aligned'])),abs(old['HedgesG_SE']-float(r['HedgesG_SE'])))
 checks.append(dict(PMID=r['PMID'],Gene=r['GeneSymbol'],OldG=r['HedgesG_Aligned'],FrozenAliasReproducedG=old['HedgesG'],ExactSymbolG=direct['HedgesG'],MaxOriginalMetricDifference=err,OriginalMechanismReproduced=err<1e-6))
assert all(r['OriginalMechanismReproduced'] for r in checks)
ps={r['FrozenLibraryAlias']:r for r in sm if r['FinalPhenotype']=='PSORIASIS_LESIONAL_SKIN'}
assert set(ps)=={f'PSO{i}' for i in range(1,8)}
assert {r['FrozenLibraryAlias'] for r in sm if r['FinalPhenotype']=='HEALTHY_CONTROL_SKIN'}=={f'HC{i}' for i in range(1,7)}
samples=[]
for c in header[1:]:
 isps=c in ps
 samples.append(dict(MatrixColumn=c,FrozenAliasDirectMatch=isps,AnalysisGroup='PSORIASIS_LESIONAL_SKIN' if isps else 'HEALTHY_CONTROL_SKIN',SpecificGSM=ps[c]['GSM'] if isps else '',Basis='DIRECT_FROZEN_ALIAS' if isps else 'SUPPLIED_RECOVERY_SET_COMPLEMENT',ExternalIndividualGSMVerified=False))
write_csv(a.out/'gene_mapping_audit.csv',diag);write_csv(a.out/'original_alias_mechanism_reproduction.csv',checks);write_csv(a.out/'sample_column_audit.csv',samples)
write_csv(a.out/'input_checksums.csv',[dict(Filename=f.name,SHA256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [nf,sf,raw,bundle]])
summary=dict(Targets=36,ExactSymbolsPresent=36,ExtraAliasMappingGenes=[r['Gene'] for r in diag if r['CorrectionNeeded']],OriginalResultsReproduced=len(checks),MaximumDifference=max(r['MaxOriginalMetricDifference'] for r in checks),PsoriasisColumnsDirectFrozenMatches=7,HealthyColumnsAssignedBySuppliedRecovery=6,IndividualHealthyGSMsVerified=False,AdditionalMetricChangesBeyondS100A9=False,Scope='Frozen-table and raw-matrix concordance. This local check does not assign individual NL-to-GSM identities. Separately retrieved official GEO records verify the 7/6 groups and public dates; see evidence/external_evidence_decisions.csv.')
(a.out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
