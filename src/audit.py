"""Step 1 - Data audit of the ZEDGER water-quality and tank exports.

What this script does, in plain language:
  1. Loads the two raw CSV exports from data/raw/.
  2. Removes readings that were deleted ("archived") in the app.
  3. Checks for duplicated days, missing days, unusual reading times and
     implausible values (outliers), per system.
  4. Rebuilds, from the tank records, how many fish and how many cleanings
     each system had on each day.
  5. Writes cleaned tables to data/clean/, a text report to reports/, and
     one figure per parameter to figures/.

Nothing is invented: every number written comes from the exports.
Run:  python src/audit.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_anonymised"  # anonymised copies (see src/anonymise.py)
CLEAN = ROOT / "data" / "clean"
FIG = ROOT / "figures"
REP = ROOT / "reports"
for p in (CLEAN, FIG, REP):
    p.mkdir(parents=True, exist_ok=True)

# The three independent water systems used in the modelling (confirmed by the facility manager).
# System 4 and System 5 (two weeks of readings only) are removed.
SYSTEMS = ["Main System", "SA Tecniplast", "SA Aquaneering"]
# Total water volume per system in litres (facility records). Main System = 3 x 175 L.
VOLUME_L = {"Main System": 525.0, "SA Tecniplast": 175.0, "SA Aquaneering": 150.0}
# First day each system's data is used. Aquaneering: data before 8 Jul ignored (facility decision).
DATA_FROM = {"Main System": "2026-04-23", "SA Tecniplast": "2026-04-23", "SA Aquaneering": "2026-07-08"}
# First day the reconstructed fish count is trusted. Before this, tanks were still being
# entered into the app (not real stocking), so the fish count is set to unknown (NaN).
FISH_RELIABLE_FROM = {"Main System": "2026-05-22", "SA Tecniplast": "2026-04-23",
                      "SA Aquaneering": "2026-07-08"}
COLORS = {"Main System": "#2a78d6", "SA Tecniplast": "#eb6834", "SA Aquaneering": "#1baf7a"}
PARAMS = {
    "temperature": ("Temperature", "°C"),
    "conductivity": ("Conductivity", "µS/cm"),
    "ph": ("pH", "pH units"),
}

# ASSUMPTION (to be confirmed): which tank rack belongs to which water system.
# "MSRack 1/2/3" = racks of the Main System. Anything else is not mapped.
def rack_to_system(rack: str | None) -> str | None:
    if not rack:
        return None
    r = rack.replace(" ", "").lower()
    if r.startswith("msrack"):
        return "Main System"
    if r == "saaquaneering":
        return "SA Aquaneering"
    if r == "satecniplast":
        return "SA Tecniplast"
    return None


KNOWN_RACKS = ["MSRack 1", "MSRack 2", "MSRack2", "MSRack3", "MSRack 3", "SA Aquaneering",
               "SA Tecniplast", "System 6", "Microplastic System 1", "Glass Stand Alone"]


def location_rack(loc: str | None) -> str | None:
    """'MSRack 1 F-1' -> 'MSRack 1'. Unknown places (incubators etc.) -> the text before the last token."""
    if not loc:
        return None
    for k in sorted(KNOWN_RACKS, key=len, reverse=True):
        if loc.startswith(k):
            return k
    return loc.rsplit(" ", 1)[0] if " " in loc else loc


# ---------------------------------------------------------------- water quality
def load_water() -> tuple[pd.DataFrame, dict]:
    w = pd.read_csv(RAW / "water_quality_rows.csv")
    info = {"rows_raw": len(w), "archived": int(w["archived"].sum())}
    w = w[~w["archived"]].copy()
    w["date"] = pd.to_datetime(w["date"])
    info["removed_other_systems"] = int((~w["system"].isin(SYSTEMS)).sum())
    w = w[w["system"].isin(SYSTEMS)]
    start = w["system"].map(DATA_FROM).pipe(pd.to_datetime)
    info["removed_before_start"] = int((w["date"] < start).sum())
    w = w[w["date"] >= start]
    w["hour"] = w["time"].str.slice(0, 2).astype(int) + w["time"].str.slice(3, 5).astype(int) / 60

    # Duplicated days: keep the LAST entry of the day (later entry = correction).
    dup = w[w.duplicated(["system", "date"], keep=False)].sort_values(["system", "date", "created_at"])
    info["duplicates"] = dup[["system", "date", "time", "temperature", "conductivity", "ph"]]
    w = w.sort_values("created_at").drop_duplicates(["system", "date"], keep="last")
    return w.sort_values(["system", "date"]), info


# Plausible ranges for a zebrafish recirculating system (outside = almost surely an error)
PLAUSIBLE = {"temperature": (20, 32), "conductivity": (150, 2500), "ph": (6.0, 9.0)}
# Minimum size of an isolated spike to be flagged
SPIKE = {"temperature": 1.5, "conductivity": 250, "ph": 0.6}


def flag_outliers(s: pd.Series, p: str) -> pd.Series:
    """Two simple, explainable rules. We only FLAG; nothing is deleted.
    1. Out of plausible range for a zebrafish system.
    2. Isolated spike: the reading jumps away from BOTH neighbours by more than
       SPIKE[p], in the same direction, while the two neighbours agree with each
       other. A real salt addition does not look like this (it stays high and
       decays slowly), but a typing or probe error does."""
    lo, hi = PLAUSIBLE[p]
    out = (s < lo) | (s > hi)
    prev, nxt = s.shift(1), s.shift(-1)
    up = (s - prev > SPIKE[p]) & (s - nxt > SPIKE[p])
    down = (prev - s > SPIKE[p]) & (nxt - s > SPIKE[p])
    agree = (prev - nxt).abs() < SPIKE[p] / 2
    return out | ((up | down) & agree)


# ---------------------------------------------------------------- tanks
def load_tanks() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    t = pd.read_csv(RAW / "tanks_rows.csv")
    end_default = pd.Timestamp("2026-11-10")
    rows, cleans, feeds = [], [], []
    unmapped_moves = 0
    for _, rec in t.iterrows():
        d = json.loads(rec["data"])
        hist = sorted(d.get("history", []), key=lambda h: (h.get("date") or "", h.get("at") or ""))
        start = pd.to_datetime(d.get("setup_date") or (hist[0]["date"] if hist else None))
        end = pd.to_datetime(str(rec["archived_at"])[:10]) if isinstance(rec["archived_at"], str) else end_default

        # location timeline: (date, rack)
        loc = [(start, location_rack(hist[0]["to"]) if hist and hist[0]["action"].startswith("Criado") else None)]
        for h in hist:
            if h["action"] in ("Criado", "Movido", "Transferido de Incubadora", "Criado em Incubadora"):
                loc.append((pd.to_datetime(h["date"]), location_rack(h["to"])))
        loc = [x for x in loc if x[1] is not None] or [(start, d.get("rack"))]
        # last known place should match current rack field
        if loc[-1][1] != d.get("rack") and d.get("rack"):
            loc.append((end, d.get("rack")))

        # fish count timeline: reconstruct backwards from mortality events
        morts = [h for h in hist if h["action"].startswith("Mortalidade")]
        n0 = int(morts[0]["from"]) if morts and str(morts[0]["from"]).isdigit() else int(d["fish_count"])
        count_changes = [(start, n0)] + [(pd.to_datetime(h["date"]), int(h["to"])) for h in morts if str(h["to"]).isdigit()]

        days = pd.date_range(start, end - pd.Timedelta(days=1) if end > start else start, freq="D")
        if len(days) == 0:
            continue
        loc_s = pd.Series({dt: r for dt, r in loc}).sort_index()
        loc_s = loc_s[~loc_s.index.duplicated(keep="last")]
        cnt_s = pd.Series({dt: n for dt, n in count_changes}).sort_index()
        cnt_s = cnt_s[~cnt_s.index.duplicated(keep="last")]
        rack_daily = loc_s.reindex(loc_s.index.union(days)).ffill().bfill().reindex(days)
        cnt_daily = cnt_s.reindex(cnt_s.index.union(days)).ffill().bfill().reindex(days)
        for dt in days:
            rows.append({"date": dt, "tank": d["id"], "rack": rack_daily[dt],
                         "system": rack_to_system(rack_daily[dt]), "fish": cnt_daily[dt],
                         "volume_l": d.get("tank_volume")})
        for c in d.get("cleanings", []):
            r = rack_daily.get(pd.to_datetime(c["date"]), d.get("rack"))
            cleans.append({"date": pd.to_datetime(c["date"]), "tank": d["id"], "system": rack_to_system(r)})
        for f in d.get("feed_log", []) or []:
            r = rack_daily.get(pd.to_datetime(f["date"]), d.get("rack"))
            feeds.append({"date": pd.to_datetime(f["date"]), "tank": d["id"], "system": rack_to_system(r),
                          "type": f.get("type"), "size": f.get("size")})
    daily = pd.DataFrame(rows)
    info = {"tanks": len(t), "tanks_archived": int(t["archived"].sum())}
    return daily, pd.DataFrame(cleans), pd.DataFrame(feeds), info


# ---------------------------------------------------------------- main
def main() -> None:
    w, winfo = load_water()
    daily_tanks, cleans, feeds, tinfo = load_tanks()
    lines: list[str] = []
    P = lines.append

    P("# Step 1 - Data audit report\n")
    P(f"Water-quality export: {winfo['rows_raw']} rows; {winfo['archived']} were deleted in the app "
      f"(archived) and removed; {winfo['removed_other_systems']} readings from System 4/5 removed; "
      f"{winfo['removed_before_start']} Aquaneering readings before 8 Jul 2026 removed.")
    P(f"Duplicated days (same system, same date) after that: {len(winfo['duplicates'])//2} pairs; "
      "the later entry of each pair was kept.\n")
    P(winfo["duplicates"].to_string(index=False) + "\n")

    P("## Coverage per system\n")
    cov = []
    for sysname, g in w.groupby("system"):
        span = pd.date_range(g.date.min(), g.date.max(), freq="D")
        missing = span.difference(g.date)
        gaps = np.diff(g.date.values).astype("timedelta64[D]").astype(int) - 1
        wk = g.date.dt.dayofweek
        cov.append({"system": sysname, "first": g.date.min().date(), "last": g.date.max().date(),
                    "calendar_days": len(span), "days_with_reading": len(g),
                    "coverage_%": round(100 * len(g) / len(span), 1), "missing_days": len(missing),
                    "longest_gap_days": int(gaps.max()) if len(gaps) else 0,
                    "weekend_readings": int((wk >= 5).sum())})
    cov = pd.DataFrame(cov)
    P(cov.to_string(index=False) + "\n")

    P("## Reading time of day\n")
    P(w.groupby("system")["hour"].describe()[["count", "min", "25%", "50%", "75%", "max"]].round(2).to_string())
    P(f"\nShare of readings between 08:30 and 09:30: {100*((w.hour>=8.5)&(w.hour<=9.5)).mean():.1f}%; "
      f"after 12:00: {100*(w.hour>=12).mean():.1f}%\n")

    P("## Value ranges and units\n")
    P(w.groupby("system")[list(PARAMS)].agg(["min", "median", "max"]).round(2).to_string())
    P("\nUnits inferred from the ranges: temperature in °C (26-30), conductivity in µS/cm "
      "(560-1570, typical for zebrafish systems), pH in pH units (6.7-9.0). No unit column exists "
      "in the export, so this is an inference.\n")

    P("## Flagged readings (possible outliers or real jumps)\n")
    flags = []
    for sysname in SYSTEMS:
        g = w[w.system == sysname].set_index("date")
        for p in PARAMS:
            f = flag_outliers(g[p], p)
            for dt in g.index[f]:
                flags.append({"system": sysname, "date": dt.date(), "parameter": p, "value": g.loc[dt, p]})
    flags = pd.DataFrame(flags)
    P(flags.to_string(index=False) if len(flags) else "None.")
    P("")

    # day-to-day jumps in conductivity (candidate salt additions / water changes)
    P("## Largest day-to-day conductivity changes (between consecutive readings)\n")
    jumps = []
    for sysname in SYSTEMS:
        g = w[w.system == sysname].sort_values("date")
        dc = g.conductivity.diff()
        dd = g.date.diff().dt.days
        for i in dc.abs().nlargest(5).index:
            jumps.append({"system": sysname, "date": g.loc[i, "date"].date(), "change_uS": dc[i],
                          "days_since_prev": int(dd[i])})
    P(pd.DataFrame(jumps).to_string(index=False) + "\n")

    # tanks -> fish per system
    P("## Tank records\n")
    P(f"Tank export: {tinfo['tanks']} tanks ({tinfo['tanks_archived']} archived). "
      "Fish count per day was rebuilt from each tank's setup date, mortality log and archive date.")
    fish = daily_tanks.groupby(["date", "system"])["fish"].sum().unstack().reindex(columns=SYSTEMS)
    for sysname, d0 in FISH_RELIABLE_FROM.items():
        fish.loc[fish.index < pd.Timestamp(d0), sysname] = np.nan
    ntank = daily_tanks.groupby(["date", "system"])["tank"].nunique().unstack().reindex(columns=SYSTEMS)
    unm = daily_tanks[daily_tanks.system.isna()].groupby("rack")["tank"].nunique()
    P("Tanks never mapped to one of the three systems (rack name):\n" + unm.to_string())
    P("\nFish per system, monthly mean of the daily total:")
    P(fish.resample("MS").mean().round(0).to_string())
    P("\nTanks per system, monthly mean:")
    P(ntank.resample("MS").mean().round(0).to_string())
    cl = cleans.groupby([cleans.date.dt.to_period("M"), "system"]).size().unstack().reindex(columns=SYSTEMS)
    P("\nCleaning events per month (tank cleanings logged in the app):")
    P(cl.to_string())
    P("\nFeeding log: present in the app but judged unreliable by the facility manager; NOT used.")

    # ------------------------------------------------ save clean data
    keep = w[w.system.isin(SYSTEMS)][["date", "time", "system"] + list(PARAMS)]
    keep.to_csv(CLEAN / "water_quality_daily.csv", index=False)
    pd.DataFrame({"system": SYSTEMS, "volume_l": [VOLUME_L[x] for x in SYSTEMS],
                  "data_from": [DATA_FROM[x] for x in SYSTEMS],
                  "fish_reliable_from": [FISH_RELIABLE_FROM[x] for x in SYSTEMS]}).to_csv(
        CLEAN / "systems.csv", index=False)
    fish.to_csv(CLEAN / "fish_per_system_daily.csv")
    cleans.groupby(["date", "system"]).size().rename("cleanings").reset_index().to_csv(
        CLEAN / "cleanings_per_system_daily.csv", index=False)
    flags.to_csv(CLEAN / "flagged_readings.csv", index=False)
    (REP / "audit_report.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))

    # ------------------------------------------------ figures
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6})
    for p, (label, unit) in PARAMS.items():
        fig, axes = plt.subplots(3, 1, figsize=(9, 6.5), sharex=True, sharey=True)
        for ax, sysname in zip(axes, SYSTEMS):
            g = w[w.system == sysname].sort_values("date")
            full = g.set_index("date")[p].asfreq("D")  # NaN on missing days -> line breaks at gaps
            ax.plot(full.index, full.values, color=COLORS[sysname], lw=1.5)
            ax.plot(g.date, g[p], "o", ms=3, color=COLORS[sysname])
            f = flags[(flags.system == sysname) & (flags.parameter == p)] if len(flags) else []
            if len(f):
                ax.plot(pd.to_datetime(f.date), f.value, "o", ms=8, mfc="none", mec="#0b0b0b", mew=1.2,
                        label="flagged")
                ax.legend(loc="upper left", frameon=False)
            c = cleans[cleans.system == sysname].date.drop_duplicates()
            for dt in c:
                ax.axvline(dt, color="#9a9994", lw=0.5, alpha=0.5, zorder=0)
            ax.set_title(sysname, loc="left", fontsize=9, color="#52514e")
            ax.set_ylabel(unit)
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        fig.suptitle(f"{label} per system (grey lines = days with a logged tank cleaning; "
                     "breaks = missing days)", x=0.01, ha="left", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG / f"audit_{p}.png", dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 3.2))
    for sysname in SYSTEMS:
        ax.plot(fish.index, fish[sysname], color=COLORS[sysname], lw=2, label=sysname)
    ax.set_ylabel("fish (reconstructed)")
    ax.legend(frameon=False, loc="upper left")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.set_title("Fish per system over time (from tank records)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "audit_fish_count.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 2.6))
    for i, sysname in enumerate(SYSTEMS):
        g = w[w.system == sysname]
        ax.plot(g.date, [i] * len(g), "|", ms=10, color=COLORS[sysname])
    ax.set_yticks(range(3), SYSTEMS)
    ax.set_ylim(-0.6, 2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.set_title("Days with a water reading (each tick = one day)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "audit_coverage.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
