"""Create the anonymised copies of the ZEDGER exports that are published in data/raw_anonymised/.

Run ONCE, locally, on the original exports (data/raw/, never published). It removes everything
that is not needed for the analysis and could identify people or unpublished projects:
user e-mails, names, free-text notes, owner ids, genetic line / project / experiment labels,
feeding logs and reminders. What remains: tank id, rack, setup date, archive date, fish count,
tank volume, the dated history of creations/moves/mortality counts, and cleaning dates; and for
water quality: date, time, system, temperature, conductivity, pH and the app's archive flag.

Run: python src/anonymise.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW, OUT = ROOT / "data" / "raw", ROOT / "data" / "raw_anonymised"
OUT.mkdir(parents=True, exist_ok=True)

w = pd.read_csv(RAW / "water_quality_rows.csv")
w[["date", "time", "system", "temperature", "conductivity", "ph", "archived", "created_at"]].to_csv(
    OUT / "water_quality_rows.csv", index=False)

t = pd.read_csv(RAW / "tanks_rows.csv")
rows = []
for _, r in t.iterrows():
    d = json.loads(r["data"])
    keep = {
        "id": d["id"], "rack": d.get("rack"), "setup_date": d.get("setup_date"),
        "fish_count": d.get("fish_count"), "tank_volume": d.get("tank_volume"),
        "history": [{k: h.get(k) for k in ("action", "date", "at", "from", "to")} for h in d.get("history", [])],
        "cleanings": [{"date": c.get("date")} for c in d.get("cleanings", [])],
    }
    rows.append({"id": d["id"], "data": json.dumps(keep, ensure_ascii=False),
                 "archived_at": r["archived_at"], "archived": r["archived"]})
pd.DataFrame(rows).to_csv(OUT / "tanks_rows.csv", index=False)
print("written", list(OUT.iterdir()))
