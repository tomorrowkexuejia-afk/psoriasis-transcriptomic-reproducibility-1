from pathlib import Path
import subprocess,sys,argparse,json,hashlib
p=argparse.ArgumentParser();p.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();r=a.release.resolve();code=Path(__file__).resolve().parent
for entry in json.loads((r/'input_manifest.json').read_text()):
 path=r/entry['Path']
 if hashlib.file_digest(path.open('rb'),'sha256').hexdigest()!=entry['SHA256']:raise RuntimeError(f'Input identity mismatch: {path}')
for cmd in [[sys.executable,str(code/'run_corrected_analysis.py'),'--project',str(r/'inputs'),'--out',str(r/'results')],[sys.executable,str(code/'audit_mapping.py'),'--project',str(r/'inputs'),'--out',str(r/'results/mapping')],[sys.executable,str(code/'build_corrected_outputs.py'),'--release',str(r)],[sys.executable,str(code/'verify_release.py'),'--release',str(r)]]:
 subprocess.run(cmd,check=True)
