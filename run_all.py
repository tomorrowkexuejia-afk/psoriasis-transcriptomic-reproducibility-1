from pathlib import Path
import subprocess, sys
root = Path(__file__).resolve().parent / "core_analysis"
subprocess.run([sys.executable, str(root / "code" / "run_pipeline.py"), "--release", str(root)], check=True)
subprocess.run([sys.executable, str(root / "code" / "audit_impact_analysis.py")], cwd=root, check=True)
print("Core analysis and paired audit-impact analysis completed successfully.")
