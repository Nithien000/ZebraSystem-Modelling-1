"""Step 4 - Parameterisation: calibrate on April-August, validate on September-October,
and quantify parameter uncertainty.

WHAT THIS SCRIPT DOES (plain language)
1. CALIBRATION. Fits the conductivity model (Step 2) and the pH model (Step 3) using only
   readings up to 31 Aug 2026. Nothing from September onwards is seen.
2. VALIDATION = honest forecast test on 1 Sep - 31 Oct 2026. For every validation reading,
   start from the PREVIOUS measured value and let the calibrated model run forward to the day
   of this reading. The model does not know about future salt additions or other jumps, so
   those show up as errors, exactly as they would in real use. We compare with the simplest
   possible forecast, "same as the last reading" (persistence), which any model must beat to
   be useful. Inputs that are logged or measured (temperature, fish count, cleaning days,
   dosing periods) are treated as known.
3. UNCERTAINTY. Moving-block bootstrap: take the calibration residuals (model minus data),
   resample them in blocks of 5 consecutive readings (to keep their day-to-day correlation),
   add them back to the fitted curve to create 300 alternative "possible" data sets, refit
   the continuous parameters on each, and report the 2.5-97.5 % range = 95 % confidence
   interval. The inferred jump dates are kept fixed in the bootstrap, so the intervals are
   CONDITIONAL on those dates (stated as a limitation).
4. STABILITY. Refit the same model on September-October only and check whether the key
   rate parameter falls inside the calibration confidence interval.

Run: python src/calibrate_validate.py
"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import conductivity_model as cm
import ph_model as pm

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CLEAN, FIG, REP = ROOT / "data" / "clean", ROOT / "figures", ROOT / "reports"
CAL_END, VAL_START, VAL_END = "2026-08-31", "2026-09-01", "2026-10-31"
SYSTEMS = cm.SYSTEMS
COLORS = cm.COLORS
N_BOOT, BLOCK = 300, 5
R_FINE = np.round(np.geomspace(0.005, 1.0, 90), 4)
RNG = np.random.default_rng(20261110)


# ------------------------------------------------------------------ fitting helpers
def fit_conductivity(name, end=None, start=None):
    m = cm.ConductivityModel(name, end, start).fit()
    mt = cm.ConductivityModel(name, end, start).fit(use_temp=True)
    r_edge = mt.params["r"] in (cm.R_GRID[0], cm.R_GRID[-1])
    if mt.best["bic"] - m.best["bic"] < -2 and not r_edge:
        m = mt
    m.kind = "cond"
    return m


def fit_ph(name, end=None, start=None):
    m = pm.PHModel(name, end, start)
    best, terms, _, _ = pm.select_terms(m)
    m.set_fit(best)
    m.kind = "ph"
    m.terms = terms
    return m


def design(m, r, k):
    """Building blocks for fixed r, k with the jump dates fixed at the fitted ones."""
    saved = m.cand
    m.cand = np.asarray(m.jump_days, dtype=int)
    if m.kind == "cond":
        names, B, J = m.blocks(r, k, m.best["use_temp"])
    else:
        names, B, J = m.blocks(r, k, m.best["terms"])
    m.cand = saved
    return names, np.column_stack([B, J]) if J.size else B


def refit_continuous(m, y):
    """Least squares for all continuous parameters with jump dates fixed (grid over r, k)."""
    best = None
    ks = [m.best["k"]] if m.best["k"] == 0 else sorted({*cm.K_GRID, *pm.K_GRID})
    for r in R_FINE:
        for k in ks:
            names, X = design(m, float(r), float(k))
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            ok = cm.physical_ok(names, beta[: len(names)]) if m.kind == "cond" else pm.physical_ok(names, beta[: len(names)])
            if not ok:
                continue
            rss = float(np.sum((X @ beta - y) ** 2))
            if best is None or rss < best[0]:
                best = (rss, r, k, names, beta)
    rss, r, k, names, beta = best
    out = dict(zip(names, beta[: len(names)]))
    out.update(r=r, k=k)
    return out, rss


def derived(m, p):
    """Interpretable quantities we want intervals for."""
    d = {"r": p["r"], "half_life_days": np.log(2) / p["r"]}
    if m.kind == "cond":
        d["Cstar"] = p["Cstar"]
        if "Cset" in p and m.has_ctrl:
            d["Cset"] = p["Cset"]
        if "g" in p:
            d["g_temp"] = p["g"]
    else:
        d["pH_eq_or_in"] = p["pH_in"]
        if "alpha" in p:
            d["alpha"] = p["alpha"]
            d["alpha_over_r"] = p["alpha"] / p["r"]
        if "delta" in p:
            d["delta_cleaning"] = p["delta"]
        if "pH_set" in p:
            d["pH_set"] = p["pH_set"]
    return d


def bootstrap(m):
    yhat = (m.simulate_daily() if m.kind == "cond" else m.simulate())[m.obs]
    res = m.y - yhat
    n = len(res)
    y0 = m.y.copy()
    draws = []
    for _ in range(N_BOOT):
        idx = np.concatenate([np.arange(s, s + BLOCK) for s in RNG.integers(0, n - BLOCK + 1, n // BLOCK + 1)])[:n]
        ystar = yhat + res[idx]
        p, _ = refit_continuous(m, ystar)
        draws.append(derived(m, p))
    m.y = y0
    return pd.DataFrame(draws)


# ------------------------------------------------------------------ forecasting
def daily_inputs(m_full, m_cal, p):
    """lam_i, forcing_i, known impulses_i on the FULL day grid, using calibrated parameters."""
    if m_cal.kind == "cond":
        act = m_full.active if m_cal.has_ctrl else np.zeros(m_full.nd)
        tdev = m_full.temp - m_cal.temp.mean()  # g was estimated relative to the calibration mean
        lam = p["r"] + p["k"] * act
        forcing = p["r"] * (p["Cstar"] + p.get("g", 0.0) * tdev) + p["k"] * act * p.get("Cset", 0.0)
        imp = np.zeros(m_full.nd)
    else:
        act = m_full.active * ("pH_set" in p)
        lam = p["r"] + p["k"] * act
        forcing = p["r"] * p["pH_in"] - p.get("alpha", 0.0) * m_full.fpl + p["k"] * act * p.get("pH_set", 0.0)
        imp = p.get("delta", 0.0) * m_full.clean
    return lam, forcing, imp


def one_step_forecasts(m_full, m_cal):
    lam, forcing, imp = daily_inputs(m_full, m_cal, m_cal.params if m_cal.kind == "cond" else m_cal.p)
    e, g = np.exp(-lam), (1 - np.exp(-lam)) / lam
    rows = []
    for j in range(1, len(m_full.obs)):
        a, b = m_full.obs[j - 1], m_full.obs[j]
        date = m_full.days[b]
        if not (pd.Timestamp(VAL_START) <= date <= pd.Timestamp(VAL_END)):
            continue
        c = m_full.y[j - 1]
        for i in range(a, b):
            c = c * e[i] + forcing[i] * g[i]
            c += imp[i + 1]
        rows.append({"date": date, "obs": m_full.y[j], "model": c, "persistence": m_full.y[j - 1],
                     "gap_days": int(b - a)})
    return pd.DataFrame(rows)


def metrics(y, yhat):
    e = np.asarray(yhat) - np.asarray(y)
    return {"RMSE": float(np.sqrt(np.mean(e ** 2))), "MAE": float(np.mean(np.abs(e))),
            "R2": float(1 - np.sum(e ** 2) / np.sum((np.asarray(y) - np.mean(y)) ** 2))}


# ------------------------------------------------------------------ main
def main():
    out = ["# Step 4 - Calibration (Apr-Aug), validation (Sep-Oct) and parameter uncertainty\n",
           f"Calibration: readings up to {CAL_END}. Validation: {VAL_START} to {VAL_END}. "
           f"Bootstrap: {N_BOOT} moving-block resamples (block = {BLOCK} readings); intervals are 95 % "
           "percentile intervals, conditional on the inferred jump dates.\n"]
    ci_rows, val_rows, stab_rows = [], [], []
    fc_store = {}
    for var, fitter, unit in [("conductivity", fit_conductivity, "µS/cm"), ("pH", fit_ph, "pH")]:
        out.append(f"## {var}\n")
        for name in SYSTEMS:
            m_cal = fitter(name, end=CAL_END)
            m_full = cm.ConductivityModel(name) if var == "conductivity" else pm.PHModel(name)
            if var == "pH":
                m_full.p = m_cal.p
            p_cal = m_cal.params if var == "conductivity" else m_cal.p
            yhat_cal = (m_cal.simulate_daily() if var == "conductivity" else m_cal.simulate())[m_cal.obs]
            cal_met = metrics(m_cal.y, yhat_cal)

            boots = bootstrap(m_cal)
            point = derived(m_cal, p_cal)
            fc = one_step_forecasts(m_full, m_cal)
            fc_store[(var, name)] = fc
            vm, vp = metrics(fc.obs, fc.model), metrics(fc.obs, fc.persistence)

            try:
                m_val = fitter(name, end=VAL_END, start=VAL_START)
                r_val = m_val.params["r"] if var == "conductivity" else m_val.p["r"]
            except Exception:
                r_val = np.nan
            lo_r, hi_r = boots["r"].quantile([0.025, 0.975])

            out.append(f"### {name}\n")
            out.append(f"- Calibration: {len(m_cal.y)} readings ({m_cal.days[0].date()} to {m_cal.days[-1].date()}); "
                       f"in-sample RMSE {cal_met['RMSE']:.3g} {unit}, R² {cal_met['R2']:.2f}; "
                       f"{len(m_cal.jumps)} inferred jumps"
                       + (f"; terms: {', '.join(sorted(m_cal.terms)) or 'relaxation only'}" if var == "pH" else ""))
            out.append("- Parameters (point estimate [95 % CI]):")
            for kname, v in point.items():
                lo, hi = boots[kname].quantile([0.025, 0.975])
                out.append(f"    - {kname}: {v:.4g} [{lo:.4g}, {hi:.4g}]")
                ci_rows.append({"variable": var, "system": name, "parameter": kname, "estimate": v,
                                "ci_low": lo, "ci_high": hi})
            out.append(f"- Validation, next-reading forecast ({len(fc)} forecasts): model RMSE {vm['RMSE']:.3g}, "
                       f"MAE {vm['MAE']:.3g}, R² {vm['R2']:.2f} | persistence RMSE {vp['RMSE']:.3g}, "
                       f"MAE {vp['MAE']:.3g}, R² {vp['R2']:.2f} -> "
                       + ("model better" if vm["RMSE"] < 0.97 * vp["RMSE"] else
                          "persistence better" if vm["RMSE"] > 1.03 * vp["RMSE"] else "about equal"))
            inside = lo_r <= r_val <= hi_r if np.isfinite(r_val) else None
            out.append(f"- Stability: r refitted on Sep-Oct only = {r_val:.3g} per day; calibration CI "
                       f"[{lo_r:.3g}, {hi_r:.3g}] -> " + ("inside (stable)" if inside else
                                                          "outside (parameter changed or poorly determined)" if inside is not None else "refit failed") + "\n")
            val_rows.append({"variable": var, "system": name, "n_val": len(fc),
                             **{f"model_{k}": v for k, v in vm.items()},
                             **{f"persistence_{k}": v for k, v in vp.items()},
                             "cal_RMSE": cal_met["RMSE"], "cal_R2": cal_met["R2"]})
            stab_rows.append({"variable": var, "system": name, "r_cal": point["r"], "r_cal_lo": lo_r,
                              "r_cal_hi": hi_r, "r_val_refit": r_val, "inside": inside})
            print(f"done {var} {name}", flush=True)

    pd.DataFrame(ci_rows).to_csv(CLEAN / "step4_parameter_ci.csv", index=False)
    pd.DataFrame(val_rows).to_csv(CLEAN / "step4_validation_metrics.csv", index=False)
    pd.DataFrame(stab_rows).to_csv(CLEAN / "step4_stability.csv", index=False)
    pd.concat({k: v for k, v in fc_store.items()}, names=["variable", "system", "i"]).to_csv(
        CLEAN / "step4_validation_forecasts.csv")
    (REP / "step4_calibration_validation.md").write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))

    # figure: validation forecasts
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    for col, (var, unit) in enumerate([("conductivity", "µS/cm"), ("pH", "pH")]):
        for row, name in enumerate(SYSTEMS):
            ax = axes[row, col]
            fc = fc_store[(var, name)]
            ax.plot(fc.date, fc.obs, "o", ms=3.5, color=COLORS[name], label="measured")
            ax.plot(fc.date, fc.model, "-", color="#52514e", lw=1.2, label="model forecast")
            ax.plot(fc.date, fc.persistence, ":", color="#9a9994", lw=1.2, label="same as last reading")
            vm, vp = metrics(fc.obs, fc.model), metrics(fc.obs, fc.persistence)
            ax.set_title(f"{name} · {var}: RMSE model {vm['RMSE']:.3g} vs persistence {vp['RMSE']:.3g}",
                         loc="left", fontsize=8.5, color="#52514e")
            ax.set_ylabel(unit)
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    for ax in axes[-1]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    fig.suptitle("Validation (Sep-Oct 2026): next-reading forecasts from models calibrated on Apr-Aug",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "step4_validation.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
