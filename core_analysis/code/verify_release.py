"""Targeted regression checks for actual numerical/provenance defects, no upstream claims."""
from pathlib import Path
import ast,csv,json,math,statistics,collections,sys,argparse
p=argparse.ArgumentParser();p.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1]);root=p.parse_args().release.resolve()
sys.path.insert(0,str(root/'code'))
from recompute import read_csv,corrected_meta
report={}
rows=read_csv(root/'results/development_recomputed.csv')
assert len(rows)==202
assert max(float(r['MaxAbsoluteDifference']) for r in rows)<1e-6
report['Original202ArithmeticMaximumDifference']=max(float(r['MaxAbsoluteDifference']) for r in rows)
# Every current target is represented by an exact-symbol row; MIF must not enter S100A9.
nom=read_csv(root/'inputs/upload/GSE295540_target_nomenclature_FINAL_FROZEN.csv')
exact={r['NCBICurrentSymbol'].upper():r['NCBIGeneID'] for r in nom}
assert len(exact)==36 and 'S100A9' in exact and 'MIF' not in exact
hold=read_csv(root/'tables/S3_corrected_holdout_39_claims.csv')
assert len(hold)==39
s9=[r for r in hold if r['Gene']=='S100A9']
assert len(s9)==2 and all(abs(float(r['HedgesG'])-4.1133963707)<1e-6 for r in s9)
report['S100A9CorrectedRows']=len(s9)
# Confirmed source uses cannot survive in either amended synthesis.
for mode,claims,contributions in [('confirmed_use_corrections',34,180),('source_adjudicated_subset',33,176)]:
 rs=read_csv(root/'tables'/f'{mode}_cohort_inputs.csv');meta=read_csv(root/'tables'/f'{mode}_development_meta.csv')
 assert len(rs)==contributions and len(meta)==claims
 for r in rs:
  assert not(r['PMID']=='39480805' and r['GSE'] in {'GSE54456','GSE66511','GSE121212'})
  assert not(r['PMID']=='41224720' and r['GSE']=='GSE121212')
 if mode=='source_adjudicated_subset':assert not any(r['PMID'] in {'39735895'} for r in rs)
 report[mode]={'Claims':claims,'Contributions':contributions}
# Extract only the numeric functions of legacy scripts; do not execute their file-writing workflow.
inputs=read_csv(root/'results/synthesis_per_cohort_inputs_corrected.csv');groups=collections.defaultdict(list)
for r in inputs:groups[(r['PMID'],r['Gene'])].append(r)
patches=[]
for f in sorted((root/'legacy_patches').rglob('*.py')):
 src=f.read_text();tree=ast.parse(src)
 funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'reml_components','golden','golden_minimize','tau_reml','estimate_tau2_reml'}]
 if not funcs:continue
 env={'math':math,'statistics':statistics};exec(compile(ast.Module(body=funcs,type_ignores=[]),str(f),'exec'),env)
 tau=env.get('tau_reml',env.get('estimate_tau2_reml'));largest=0
 for rs in groups.values():
  y=[float(r['HedgesG_Aligned']) for r in rs];v=[float(r['HedgesG_SE'])**2 for r in rs]
  oldfn=tau(y,v);newfn=corrected_meta(rs)['Tau2'];largest=max(largest,abs(oldfn-newfn))
 assert largest<1e-4,(f,largest)
 patches.append({'File':str(f.relative_to(root)),'ComparedClaims':len(groups),'MaxTau2Difference':largest})
report['LegacyNumericalFunctions']=patches
# PowerShell embeds exactly the repaired Python implementation.
wrappers=[]
for ps in sorted((root/'legacy_patches').rglob('*.ps1')):
 txt=ps.read_text(encoding='utf-8-sig');start=txt.index("@'\n")+4;end=txt.index("\n'@",start);embedded=txt[start:end]
 matches=[f for f in ps.parent.glob('*.py') if "'"+f.name+"'" in txt]
 assert len(matches)==1
 assert ast.dump(ast.parse(embedded))==ast.dump(ast.parse(matches[0].read_text()))
 wrappers.append(ps.name)
report['PowerShellEmbeddedPythonMatches']=wrappers
report['Scope']='Portable analysis outputs and extracted legacy numeric functions verified. Full legacy staged integration is not run because its upstream artifacts are incomplete.'
(root/'evidence/regression_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
