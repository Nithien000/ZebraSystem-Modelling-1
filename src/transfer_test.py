"""Step 5 - Transfer test: does the model learned on Main System and SA Tecniplast work on
the independent system SA Aquaneering (running since July 2026)?

PLAIN-LANGUAGE DESIGN
For each donor system (Main, Tecniplast) and each variable (conductivity, pH) we compare
three ways of modelling Aquaneering:
  (A) FULL TRANSFER  - take the donor's process parameters unchanged (renewal/relaxation rate r,
      equilibrium C* or pH_in, temperature or cleaning effect if the donor had one).
  (B) RATE TRANSFER  - take only the donor's rate r; estimate Aquaneering's own equilibrium.
  (C) REFIT          - estimate everything on Aquaneering itself (Steps 2-4).
Events that are specific to a system (inferred salt additions / dilutions / other jumps) are
always re-inferred for Aquaneering, because a donor cannot know what was done to another system.

Two tests:
  1. Fit to the whole Aquaneering record (8 Jul - 10 Nov): RMSE and BIC. BIC only counts the
     parameters actually estimated on Aquaneering, so transferred parameters are "free".
  2. Forecast test on Sep-Oct (as in Step 4): next-reading forecasts. For a clean time split,
     the donor parameters here come from the donors' Apr-Aug calibrations, and Aquaneering's own
     parameters (B, C) from its Jul-Aug calibration.

Run: python src/transfer_test.py
"""
from __future__ import annotations

import warnings
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import calibrate_validate as cv
import conductivity_model as cm
import ph_model as pm

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CLEAN, FIG, REP = ROOT / "data" / "clean", ROOT / "figures", ROOT / "reports"
TARGET = "SA Aquaneering"
DONORS = ["Main System", "SA Tecniplast"]
COLORS = cm.COLORS


def params_of(m):
    return m.params if m.kind == "cond" else m.p


def transferable(m, rate_only: bool):
    """Which donor parameters are carried over. Dosing set-points are NOT transferred:
    Aquaneering shows no dosing periods."""
    p = params_of(m)
    if m.kind == "cond":
        keys = ["r"] if rate_only else ["r", "Cstar"] + (["g"] if "g" in p else [])
    else:
        keys = ["r"] if rate_only else ["r", "pH_in"] + (["delta"] if "delta" in p else [])
    return {k: float(p[k]) for k in keys}


def fit_fixed(m, fixed: dict):
    """Fit model m with some parameters FIXED (r always fixed). Jumps are re-inferred
    (LASSO + modified BIC, as in Steps 2-3). Returns predictions, rss, bic, n_jumps, params."""
    r = fixed["r"]
    if m.kind == "cond":
        arg = "g" in fixed
    else:
        arg = {"cleaning"} if "delta" in fixed else set()
    orig = m.blocks
    has_ctrl0 = getattr(m, "has_ctrl", False)
    if m.kind == "cond":
        m.has_ctrl = False  # set before building the column list, so the kept indices match
    names0, B0, _ = orig(r, 0.0, arg)
    off = sum(B0[:, i] * fixed[n] for i, n in enumerate(names0) if n in fixed) if any(
        n in fixed for n in names0) else np.zeros(len(m.y))
    keep = [i for i, n in enumerate(names0) if n not in fixed]

    def blocks(r_, k_, a_):
        names, B, J = orig(r_, k_, a_)
        return [names[i] for i in keep], B[:, keep], J

    y0 = m.y.copy()
    m.blocks, m.y = blocks, y0 - off
    best = m.fit_given(r, 0.0, arg)
    m.blocks, m.y = orig, y0
    names, B, J = blocks(r, 0.0, arg)
    if m.kind == "cond":
        m.has_ctrl = has_ctrl0
    sup = list(best["sup"])
    X = np.column_stack([B, J[:, sup]]) if sup else B
    pred = X @ best["beta"] + off
    rss = float(np.sum((pred - y0) ** 2))
    n = len(y0)
    nfree = len(keep)
    bic = n * np.log(rss / n) + (nfree + cm.JUMP_PENALTY * len(sup)) * np.log(n)
    p = dict(fixed)
    p.update(dict(zip(names, best["beta"][: len(names)])))
    p.setdefault("k", 0.0)
    return pred, rss, bic, len(sup), p


def refit_bic(m):
    """BIC of the own-fit (C), counted the same way (estimated params incl. r)."""
    pred = (m.simulate_daily() if m.kind == "cond" else m.simulate())[m.obs]
    rss = float(np.sum((pred - m.y) ** 2))
    n = len(m.y)
    nb = len(m.best["names"]) + 1
    return pred, rss, n * np.log(rss / n) + (nb + cm.JUMP_PENALTY * len(m.jumps)) * np.log(n), len(m.jumps)


def fake_cal(kind, p, temp_ref):
    """Minimal object so Step 4's forecast function can use any parameter set."""
    p = {**p}
    p.setdefault("k", 0.0)
    ns = SimpleNamespace(kind=kind, has_ctrl=False, temp=temp_ref)
    if kind == "cond":
        ns.params = p
    else:
        ns.p = p
    return ns


def main():
    out = ["# Step 5 - Transfer test on SA Aquaneering\n",
           "(A) full transfer of donor process parameters, (B) donor rate only, (C) Aquaneering's own fit.",
           "System-specific jumps are always re-inferred on Aquaneering.\n"]
    rows = []
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.5), sharex=True)
    for ax, (var, fitter, unit) in zip(axes, [("conductivity", cv.fit_conductivity, "µS/cm"),
                                             ("pH", cv.fit_ph, "pH")]):
        out.append(f"## {var}\n")
        own_full = fitter(TARGET)
        pred_c, rss_c, bic_c, nj_c = refit_bic(own_full)
        own_cal = fitter(TARGET, end=cv.CAL_END)
        m_full_fc = cm.ConductivityModel(TARGET) if var == "conductivity" else pm.PHModel(TARGET)
        temp_ref = own_cal.temp if var == "conductivity" else None
        fc_c = cv.one_step_forecasts(m_full_fc, own_cal)
        p_own = params_of(own_full)
        out.append(f"- (C) Aquaneering own fit: r = {p_own['r']:.3g}, "
                   + (f"C* = {p_own['Cstar']:.0f}" if var == "conductivity" else f"pH_eq = {p_own['pH_in']:.2f}")
                   + f"; RMSE {np.sqrt(rss_c/len(own_full.y)):.3g} {unit}, BIC {bic_c:.1f}, {nj_c} jumps; "
                   f"forecast RMSE {cv.metrics(fc_c.obs, fc_c.model)['RMSE']:.3g} "
                   f"(persistence {cv.metrics(fc_c.obs, fc_c.persistence)['RMSE']:.3g})")
        rows.append({"variable": var, "donor": "-", "variant": "C refit", "r": p_own["r"],
                     "equilibrium": p_own.get("Cstar", p_own.get("pH_in")),
                     "fit_RMSE": np.sqrt(rss_c / len(own_full.y)), "BIC": bic_c, "dBIC_vs_C": 0.0, "jumps": nj_c,
                     "forecast_RMSE": cv.metrics(fc_c.obs, fc_c.model)["RMSE"],
                     "persistence_RMSE": cv.metrics(fc_c.obs, fc_c.persistence)["RMSE"]})
        ax.plot(own_full.g.date, own_full.y, "o", ms=3.5, color=COLORS[TARGET], label="measured")
        ax.plot(own_full.g.date, pred_c, "-", color="#52514e", lw=1.4, label="(C) own fit")

        for donor in DONORS:
            d_full, d_cal = fitter(donor), fitter(donor, end=cv.CAL_END)
            for variant, rate_only in [("A full transfer", False), ("B rate only", True)]:
                # 1) whole-record fit with donor (full-period) parameters
                m = cm.ConductivityModel(TARGET) if var == "conductivity" else pm.PHModel(TARGET)
                m.kind = "cond" if var == "conductivity" else "ph"
                fixed = transferable(d_full, rate_only)
                pred, rss, bic, nj, p = fit_fixed(m, fixed)
                # 2) forecast test with donor CALIBRATION parameters (clean time split)
                fixed_cal = transferable(d_cal, rate_only)
                mc = cm.ConductivityModel(TARGET, end=cv.CAL_END) if var == "conductivity" else pm.PHModel(TARGET, end=cv.CAL_END)
                mc.kind = m.kind
                _, _, _, _, p_cal = fit_fixed(mc, fixed_cal)
                if not rate_only:
                    p_cal = {**p_cal, **fixed_cal}
                fc = cv.one_step_forecasts(m_full_fc, fake_cal(m.kind, p_cal, mc.temp if var == "conductivity" else None))
                fm = cv.metrics(fc.obs, fc.model)["RMSE"]
                eq = p.get("Cstar", p.get("pH_in"))
                out.append(f"- ({variant[0]}) from {donor}: r = {p['r']:.3g} "
                           + "(transferred)"
                           + f", equilibrium = {eq:.4g} " + ("(transferred)" if not rate_only else "(Aquaneering's own)")
                           + f"; RMSE {np.sqrt(rss/len(m.y)):.3g} {unit}, BIC {bic:.1f} (Δ vs C {bic-bic_c:+.1f}), "
                           f"{nj} jumps; forecast RMSE {fm:.3g}")
                rows.append({"variable": var, "donor": donor, "variant": variant, "r": p["r"], "equilibrium": eq,
                             "fit_RMSE": np.sqrt(rss / len(m.y)), "BIC": bic, "dBIC_vs_C": bic - bic_c, "jumps": nj,
                             "forecast_RMSE": fm, "persistence_RMSE": cv.metrics(fc.obs, fc.persistence)["RMSE"]})
                if not rate_only:
                    ax.plot(m.g.date, pred, "--", lw=1.1, color=COLORS[donor],
                            label=f"(A) parameters from {donor}")
            print(f"done {var} {donor}", flush=True)
        ax.set_ylabel(unit)
        ax.set_title(f"SA Aquaneering {var}: own fit vs full transfer from the other systems",
                     loc="left", fontsize=9, color="#52514e")
        if var == "conductivity":
            ax.legend(frameon=False, fontsize=7.5, loc="upper right")
        out.append("")
    import matplotlib.dates as mdates
    axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.tight_layout()
    fig.savefig(FIG / "step5_transfer.png", dpi=150)
    plt.close(fig)
    tab = pd.DataFrame(rows)
    tab.to_csv(CLEAN / "step5_transfer.csv", index=False)
    out += ["## Summary table\n", tab.round(3).to_string(index=False)]
    (REP / "step5_transfer.md").write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
