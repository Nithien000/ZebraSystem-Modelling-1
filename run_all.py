"""Reproduce every result, table and figure of the project, in order.

    pip install -r requirements.txt
    python run_all.py

Takes about 5 minutes on a laptop. Inputs: data/raw_anonymised/. Outputs: data/clean/,
figures/, reports/.
"""
import runpy
import sys
import time
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
sys.path.insert(0, str(SRC))
STEPS = [
    ("1 data audit", "audit.py"),
    ("2 conductivity model", "conductivity_model.py"),
    ("3 pH model", "ph_model.py"),
    ("4 calibration / validation / bootstrap", "calibrate_validate.py"),
    ("5 transfer test", "transfer_test.py"),
    ("6 ML comparison", "ml_compare.py"),
]
for label, script in STEPS:
    t0 = time.time()
    print(f"\n=== Step {label} ===", flush=True)
    runpy.run_path(str(SRC / script), run_name="__main__")
    print(f"=== Step {label} done in {time.time() - t0:.0f} s ===", flush=True)
