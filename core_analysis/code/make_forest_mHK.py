"""Claim-level forest plot of the source-adjudicated 33-claim subset drawn with the
modified Hartung-Knapp 95% confidence interval (the manuscript's primary inference).

build_corrected_outputs.py draws the same figure with the conventional normal 95% CI
(FigureS1_Adjudicated_Forest.*).  This companion script uses identical data, identical
ordering and identical plotting parameters; only the interval columns differ
(CI_Lower/CI_Upper -> mHK_Lower/mHK_Upper) and the axis label names the interval type.
It produces the manuscript's Figure S1.

Usage:  python make_forest_mHK.py --release <path to core_analysis>
"""
from pathlib import Path
import argparse, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);root=p.parse_args().release.resolve()
sys.path.insert(0,str(root/'code'))
from recompute import read_csv
fig=root/'figures'
if not fig.is_dir():fig.mkdir(parents=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(f,n):
 for ext in ['png','pdf','svg']:f.savefig(fig/(n+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(f)

r=sorted(read_csv(root/'tables'/'source_adjudicated_subset_development_meta.csv'),key=lambda x:float(x['HedgesG']))
assert len(r)==33,len(r)
above=sum(1 for x in r if float(x['mHK_Lower'])>0)
below=sum(1 for x in r if float(x['mHK_Upper'])<0)
cross=[x['Gene'] for x in r if float(x['mHK_Lower'])<0<float(x['mHK_Upper'])]
# Section 3.5 of the manuscript reports exactly these classifications.
assert above==30 and below==0 and len(cross)==3,(above,below,cross)
assert sorted(cross)==sorted(['CD28','ITGAL','TLN1']),cross

f,ax=plt.subplots(figsize=(8,11))
for i,x in enumerate(r):
 g=float(x['HedgesG']);ax.errorbar(g,i,xerr=[[g-float(x['mHK_Lower'])],[float(x['mHK_Upper'])-g]],fmt='o',color='#356e8d',ms=4,capsize=2)
ax.set_yticks(range(len(r)),[f"{x['Gene']} | PMID {x['PMID']}" for x in r],fontsize=8)
ax.axvline(0,color='black',lw=.8)
ax.set_xlabel('Pooled aligned Hedges g (modified Hartung-Knapp 95% CI)');ax.set_title('Source-adjudicated subset: 33 claims',loc='left')
save(f,'FigureS1_Adjudicated_Forest_mHK')
print(f"claims={len(r)} mHK entirely above zero={above} entirely below zero={below} crossing zero={cross}")
