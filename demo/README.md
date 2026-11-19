# Demo on synthetic data

> **Disclaimer.** Everything in this folder is made from **invented (synthetic) data**. The real
> project used daily water-quality records from a research zebrafish facility. Those data belong to
> the host institution and **cannot be made public**, so they are not in this repository. Figures
> and reports in `demo/workspace/` say nothing about any real system.

## What the demo does

1. `make_synthetic_data.py` simulates three recirculating systems for 23 Apr to 10 Nov 2026 from
   known equations and known, invented parameters. It adds measurement noise, missing days
   (about 38 %), a few typing errors, salt additions, dilutions, fish stocking and tank cleanings,
   and writes tables in the same format as the real exports.
2. `run_demo.py` runs the **unchanged** project code (steps 1-6) on those tables and stamps every
   figure and report as synthetic.

```bash
pip install -r requirements.txt
python demo/run_demo.py      # about 4 minutes
```

## Can the fitting code recover parameters it was given?

Because the true values are known, the demo is also a test of the method (step 2, conductivity):

| | True renewal *r* (per day) | Fitted *r* | True *C\** (µS/cm) | Fitted *C\** | True jumps | Jumps found |
|---|---|---|---|---|---|---|
| Main System | 0.12 | 0.129 | 640 | 638 | 3 | 3 (dates within 1 day) |
| SA Tecniplast | 0.20 | 0.209 | 720 | 709 | 2 | 2 (within 1 day) |
| SA Aquaneering | 0.05 | 0.044 | 1000 | 986 | 1 | 1 (2 days early) |

Renewal rates, equilibria and salt-addition events are recovered well. **The automatic-dosing
detector is not**: it missed the one real dosing period (Tecniplast, May-Jul) and flagged short
flat stretches elsewhere as dosing. That detector is a simple rule (5+ readings with low
variability), and the demo shows it should be replaced by controller logs where available.
