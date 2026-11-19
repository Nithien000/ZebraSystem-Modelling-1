"""Step 6 - Machine-learning comparison: next-reading forecasts of conductivity and pH.

PLAIN-LANGUAGE DESIGN
Task: at each reading, predict the NEXT reading of the same system (usually the next day;
after a weekend or gap, a few days later). Exactly the same task and the same test points as
the mechanistic forecasts in Step 4, so the comparison is fair.

No leakage:
  * Training data = targets up to 31 Aug 2026. Test data = targets 1 Sep - 31 Oct 2026.
  * Every feature for a target uses only information available BEFORE that reading:
    previous readings, the gap length, the previous temperature, fish per litre and whether
    a cleaning was logged in the gap.
  * Model settings (ridge penalty) are chosen with time-series cross-validation inside the
    training period only (always training on the past, validating on the future).

Methods (kept small: ~200 training rows in total):
  1. Persistence          - next reading = last reading (the baseline every method must beat).
  2. Ridge regression     - linear model on lags, predicting the CHANGE from the last reading.
  3. Gradient boosting    - small tree ensemble (depth 2, 150 trees, learning rate 0.05),
                            also predicting the change.
  For 2 and 3 we try per-system models AND one pooled model over the three systems (with
  system indicators), because pooling triples the training data.
  4. Mechanistic model    - the Step 4 forecasts (calibrated ODE run from the last reading).

Uncertainty: 95 % bootstrap interval (2000 resamples of test points) of the RMSE difference
between each method and persistence. If the interval contains 0, the difference is not
distinguishable from chance with this little data.

Run: python src/ml_compare.py
"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import TimeSeriesSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CLEAN, FIG, REP = ROOT / "data" / "clean", ROOT / "figures", ROOT / "reports"
SYSTEMS = ["Main System", "SA Tecniplast", "SA Aquaneering"]
CAL_END, VAL_START, VAL_END = pd.Timestamp("2026-08-31"), pd.Timestamp("2026-09-01"), pd.Timestamp("2026-10-31")
VARS = {"conductivity": "conductivity", "pH": "ph"}
RNG = np.random.default_rng(6)
SEED = 6


def build_rows(var_col: str) -> pd.DataFrame:
    w = pd.read_csv(CLEAN / "water_quality_daily.csv", parse_dates=["date"])
    fish = pd.read_csv(CLEAN / "fish_per_system_daily.csv", parse_dates=["date"]).set_index("date")
    vol = pd.read_csv(CLEAN / "systems.csv").set_index("system").volume_l
    cl = pd.read_csv(CLEAN / "cleanings_per_system_daily.csv", parse_dates=["date"])
    # the same isolated spikes flagged in the audit are left out, as in the mechanistic models
    fl = pd.read_csv(CLEAN / "flagged_readings.csv", parse_dates=["date"])
    fl = fl[fl.parameter == var_col]
    rows = []
    for s in SYSTEMS:
        g = w[w.system == s].sort_values("date")
        g = g[~g.date.isin(fl[fl.system == s].date)].reset_index(drop=True)
        y, d, T = g[var_col].values, g.date.values, g.temperature.values
        f = fish[s].ffill()
        cls = cl[cl.system == s].set_index("date").cleanings
        for i in range(3, len(g)):
            prev, now = pd.Timestamp(d[i - 1]), pd.Timestamp(d[i])
            fpl = f.get(prev, np.nan) / vol[s]
            rows.append({
                "system": s, "date": now, "y": y[i],
                "lag1": y[i - 1], "lag2": y[i - 2], "lag3": y[i - 3],
                "d1": y[i - 1] - y[i - 2], "d2": y[i - 2] - y[i - 3],
                "gap": (now - prev).days, "gap_prev": (prev - pd.Timestamp(d[i - 2])).days,
                "temp_prev": T[i - 1], "fish_per_l": fpl,
                "cleaned_in_gap": float(cls[(cls.index > prev) & (cls.index <= now)].sum() > 0),
            })
    df = pd.DataFrame(rows)
    df["fish_per_l"] = df.groupby("system")["fish_per_l"].transform(lambda x: x.fillna(x.median()))
    for s in SYSTEMS:
        df[f"sys_{s.split()[-1]}"] = (df.system == s).astype(float)
    return df


FEATS = ["lag1", "lag2", "lag3", "d1", "d2", "gap", "gap_prev", "temp_prev", "fish_per_l", "cleaned_in_gap"]
SYS_FEATS = ["sys_System", "sys_Tecniplast", "sys_Aquaneering"]


def ridge():
    return make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 30), cv=TimeSeriesSplit(4)))


def gbr():
    return GradientBoostingRegressor(max_depth=2, n_estimators=150, learning_rate=0.05,
                                     subsample=0.8, random_state=SEED)


def fit_predict(make, train, test, feats):
    """Models predict the CHANGE from the last reading; we add it back."""
    tr = train.sort_values("date")
    m = make().fit(tr[feats].values, (tr.y - tr.lag1).values)
    return test.lag1.values + m.predict(test[feats].values)


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def boot_diff(y, pred, pers, n=2000):
    idx = RNG.integers(0, len(y), (n, len(y)))
    d = np.sqrt(((pred[idx] - y[idx]) ** 2).mean(1)) - np.sqrt(((pers[idx] - y[idx]) ** 2).mean(1))
    return np.percentile(d, [2.5, 97.5])


def main():
    mech = pd.read_csv(CLEAN / "step4_validation_forecasts.csv", parse_dates=["date"])
    out = ["# Step 6 - ML vs mechanistic vs persistence (next-reading forecasts, test Sep-Oct 2026)\n",
           "RMSE on the test period; in brackets the 95 % bootstrap interval of (method RMSE - persistence RMSE).",
           "Negative = better than persistence. An interval containing 0 = not distinguishable from persistence.\n"]
    res = []
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.2), sharey=False)
    for r_i, (var, col) in enumerate(VARS.items()):
        df = build_rows(col)
        train_all = df[df.date <= CAL_END]
        test_all = df[(df.date >= VAL_START) & (df.date <= VAL_END)]
        out.append(f"## {var}\n")
        pooled = {
            "Ridge (pooled)": fit_predict(ridge, train_all, test_all, FEATS + SYS_FEATS),
            "Gradient boosting (pooled)": fit_predict(gbr, train_all, test_all, FEATS + SYS_FEATS),
        }
        for c_i, s in enumerate(SYSTEMS):
            tr, te = train_all[train_all.system == s], test_all[test_all.system == s]
            mask = (test_all.system == s).values
            preds = {
                "Persistence": te.lag1.values,
                "Ridge (per system)": fit_predict(ridge, tr, te, FEATS),
                "Gradient boosting (per system)": fit_predict(gbr, tr, te, FEATS),
                **{k: v[mask] for k, v in pooled.items()},
            }
            mk = mech[(mech.variable == var) & (mech.system == s)].set_index("date").model
            mech_pred = mk.reindex(te.date).values
            if np.isfinite(mech_pred).all():
                preds["Mechanistic ODE (Step 4)"] = mech_pred
            y = te.y.values
            pers = preds["Persistence"]
            out.append(f"### {s}  (train {len(tr)} rows, test {len(te)} rows)\n")
            for name, p in preds.items():
                e = rmse(y, p)
                lo, hi = boot_diff(y, p, pers) if name != "Persistence" else (0.0, 0.0)
                verdict = "" if name == "Persistence" else (
                    "better" if hi < 0 else "worse" if lo > 0 else "no clear difference")
                out.append(f"- {name}: RMSE {e:.3g}" + ("" if name == "Persistence" else
                                                        f"  [{lo:+.3g}, {hi:+.3g}] -> {verdict}"))
                res.append({"variable": var, "system": s, "method": name, "RMSE": e,
                            "MAE": float(np.mean(np.abs(y - p))), "diff_lo": lo, "diff_hi": hi,
                            "verdict": verdict, "n_train": len(tr), "n_test": len(te)})
            out.append("")
            # figure: RMSE relative to persistence
            ax = axes[r_i, c_i]
            sub = pd.DataFrame([r for r in res if r["variable"] == var and r["system"] == s])
            rel = sub.RMSE / sub.loc[sub.method == "Persistence", "RMSE"].values[0]
            colors = ["#9a9994" if m == "Persistence" else "#52514e" if m.startswith("Mech") else "#2a78d6"
                      for m in sub.method]
            ax.barh(range(len(sub)), rel, color=colors, height=0.6)
            ax.axvline(1, color="#0b0b0b", lw=0.8)
            ax.set_yticks(range(len(sub)), sub.method if c_i == 0 else [""] * len(sub), fontsize=7.5)
            ax.invert_yaxis()
            ax.set_title(f"{s} · {var}", loc="left", fontsize=8.5, color="#52514e")
            ax.set_xlim(0.6, max(1.6, rel.max() + 0.05))
            if r_i == 1:
                ax.set_xlabel("RMSE / persistence RMSE  (<1 = better)", fontsize=8)
    fig.suptitle("Next-reading forecasts, test Sep-Oct 2026: error relative to 'same as last reading'",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "step6_ml_comparison.png", dpi=150)
    plt.close(fig)
    tab = pd.DataFrame(res)
    tab.to_csv(CLEAN / "step6_ml_comparison.csv", index=False)
    (REP / "step6_ml_comparison.md").write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
