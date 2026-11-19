"""Create SYNTHETIC (invented) water-quality and tank tables in the same format as the real
facility exports, so the whole pipeline can be demonstrated without any real data.

Everything here is simulated from a known set of equations with the parameters listed below.
None of the numbers are real measurements. Because the true parameters are known, the demo also
shows whether the fitting code can recover them.

    python demo/make_synthetic_data.py <output_folder>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
START, END = pd.Timestamp("2026-04-23"), pd.Timestamp("2026-11-10")

# True (invented) parameters. r = fraction of the water renewed per day.
TRUE = {
    "Main System":    dict(r=0.12, Cstar=640.0, pH_in=7.80, rph=0.10, alpha=0.004, vol=525.0, start=START),
    "SA Tecniplast":  dict(r=0.20, Cstar=720.0, pH_in=7.60, rph=0.15, alpha=0.004, vol=175.0, start=START),
    "SA Aquaneering": dict(r=0.05, Cstar=1000.0, pH_in=7.70, rph=0.08, alpha=0.0, vol=150.0,
                           start=pd.Timestamp("2026-07-08")),
}
# Invented salt additions (+) and dilutions (-) in uS/cm, by date.
JUMPS = {
    "Main System": {"2026-06-10": 250, "2026-08-02": -180, "2026-09-20": 220},
    "SA Tecniplast": {"2026-08-20": 260, "2026-10-05": -200},
    "SA Aquaneering": {"2026-09-01": 300},
}
# Invented automatic-dosing window (controller holds conductivity near a set-point).
DOSING = {"SA Tecniplast": ("2026-05-06", "2026-07-29", 700.0, 0.6)}
RACK = {"Main System": "MSRack 1", "SA Tecniplast": "SA Tecniplast", "SA Aquaneering": "SA Aquaneering"}
# (system, number of tanks, fish per tank, setup date)
TANKS = [("Main System", 8, 160, "2026-05-01"), ("SA Tecniplast", 4, 30, "2026-04-20"),
         ("SA Tecniplast", 6, 35, "2026-09-14"), ("SA Aquaneering", 10, 80, "2026-07-01")]


def main(out: Path) -> None:
    rng = np.random.default_rng(SEED)
    out.mkdir(parents=True, exist_ok=True)
    days = pd.date_range(START, END, freq="D")
    # ---- tanks
    tank_rows, fish = [], {s: np.zeros(len(days)) for s in TRUE}
    n = 1000
    for system, k, per, setup in TANKS:
        for _ in range(k):
            n += 1
            tid = f"DEMO-{n}"
            setup_d = pd.Timestamp(setup)
            cleans = [str((setup_d + pd.Timedelta(days=int(x))).date())
                      for x in range(25, (END - setup_d).days, int(rng.integers(28, 45)))]
            d = {"id": tid, "rack": RACK[system], "setup_date": setup, "fish_count": per,
                 "tank_volume": 3.5,
                 "history": [{"action": "Criado", "date": setup, "at": None, "from": None,
                              "to": f"{RACK[system]} A-{n % 9 + 1}"}],
                 "cleanings": [{"date": c} for c in cleans]}
            tank_rows.append({"id": tid, "data": json.dumps(d), "archived_at": None, "archived": False})
            fish[system][days >= setup_d] += per
    pd.DataFrame(tank_rows).to_csv(out / "tanks_rows.csv", index=False)

    # ---- water readings
    rows = []
    for system, p in TRUE.items():
        C, pH = p["Cstar"], p["pH_in"]
        for i, day in enumerate(days):
            temp_mean = 27.5 + 1.2 * np.sin((day.dayofyear - 120) / 365 * 2 * np.pi)
            jump = JUMPS.get(system, {}).get(str(day.date()), 0.0)
            dose = 0.0
            if system in DOSING:
                a, b, cset, k = DOSING[system]
                if pd.Timestamp(a) <= day <= pd.Timestamp(b):
                    dose = k * (cset - C)
            C = C + p["r"] * (p["Cstar"] - C) + dose + jump
            fpl = fish[system][i] / p["vol"]
            pH = pH + p["rph"] * (p["pH_in"] - pH) - p["alpha"] * fpl + rng.normal(0, 0.03)
            if day < p["start"] or rng.random() > 0.62:
                continue
            c_obs = C + rng.normal(0, 25)
            ph_obs = pH + rng.normal(0, 0.10)
            if rng.random() < 0.015:                      # an invented typing/probe error
                c_obs += rng.choice([-1, 1]) * 400
            hh, mm = int(rng.integers(10, 13)), int(rng.integers(0, 60))
            rows.append({"date": str(day.date()), "time": f"{hh:02d}:{mm:02d}", "system": system,
                         "temperature": round(temp_mean + rng.normal(0, 0.4), 1),
                         "conductivity": round(c_obs), "ph": round(ph_obs, 1), "archived": False,
                         "created_at": f"{day.date()} {hh:02d}:{mm:02d}:00.000+00"})
    pd.DataFrame(rows).to_csv(out / "water_quality_rows.csv", index=False)
    print(f"Wrote {len(rows)} synthetic readings and {len(tank_rows)} synthetic tanks to {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
