"""Step 2 - Process-based (mechanistic) model of conductivity, one per water system.

PLAIN-LANGUAGE IDEA
Conductivity measures how much dissolved salt is in the water. In a recirculating
system the salt amount changes for four reasons:
  * renewal: every day a fraction r of the water is drained and replaced with makeup
    water (RO water + a little salt/buffer). This pulls conductivity towards the makeup
    level;
  * evaporation: pure water leaves, salts stay; the level is topped up with makeup
    water, which adds a little more salt. Together with renewal this sets an
    equilibrium value C* that the system drifts to if nobody intervenes;
  * automatic dosing: when the controller works, it adds salt whenever conductivity
    falls below its set-point C_set, holding the system near C_set;
  * manual salt additions (not logged): sudden upward jumps. Big unlogged water
    changes look like sudden downward jumps.

EQUATION (ordinary differential equation; C in µS/cm, t in days)

    dC/dt = r * (C* + g * (T - T_mean) - C)        renewal + evaporation (g: optional temperature effect)
          + a(t) * k * (C_set - C)                 automatic dosing; a(t) = 1 when the controller is ON
          + sum_j  D_j * delta(t - t_j)            jumps: inferred salt additions (+) or dilutions (-)

Why C* and not separate "makeup conductivity" and "evaporation rate": with evaporation
E (L/day) topped up with makeup water of conductivity C_in, the equilibrium is
C* = C_in * (1 + E / (r V)). The data only show where conductivity settles (C*) and how fast
it gets there (r). C_in and E cannot be separated without measuring the makeup water.

FITTING (deterministic, explainable)
For a given renewal rate r (and dosing strength k) the model output is a straight-line
(linear) combination of a few building blocks: the starting value, C*, g, C_set and one
possible jump in every gap between two readings. So:
  1. For each r on a grid (0.005 to 1 per day) and each k on a grid, we compute those
     building blocks by solving the ODE exactly day by day.
  2. A sparse regression (LASSO) chooses which gaps really need a jump: every gap is
     allowed one, but each jump costs a penalty, so only jumps the data clearly need survive.
  3. For every candidate set of jumps we refit by ordinary least squares and compute the
     Bayesian Information Criterion (BIC = misfit + penalty for the number of parameters).
     The lowest BIC wins, across all r, k and jump sets.
  4. The final model is re-simulated with scipy.integrate.solve_ivp to check it.
We also refit with r fixed at 0.10 per day (the facility's 10 %/day renewal) and compare BIC.
Jumps are placed on the middle day of the gap in which they happened (timing ASSUMPTION).

Run: python src/conductivity_model.py
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
from scipy.integrate import solve_ivp
from sklearn.linear_model import lasso_path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CLEAN, FIG, REP = ROOT / "data" / "clean", ROOT / "figures", ROOT / "reports"
SYSTEMS = ["Main System", "SA Tecniplast", "SA Aquaneering"]
COLORS = {"Main System": "#2a78d6", "SA Tecniplast": "#eb6834", "SA Aquaneering": "#1baf7a"}
MG_PER_L_PER_US = 0.5      # rough NaCl-like conversion: 1 µS/cm ~ 0.5 mg/L (ASSUMPTION)
STABLE_SD = 15.0           # µS/cm; rolling SD below this over >= 5 readings -> controller ON
STABLE_MIN_READINGS = 5
MIN_JUMP = 75.0            # µS/cm (~10 % of typical level); smaller changes treated as noise (ASSUMPTION)
JUMP_PENALTY = 3.0         # each jump costs 3*ln(n) in the BIC: 1 for its size + 2 for choosing WHERE
                           # it goes among ~100 possible gaps (standard "modified BIC" for change points)
R_GRID = np.round(np.geomspace(0.005, 1.0, 45), 4)
K_GRID = [0.5, 1.0, 2.0, 5.0, 10.0]
# Physically plausible ranges. Fits giving values outside them are rejected (ASSUMPTION):
# an equilibrium or set-point below RO-makeup levels or above brackish levels is not realistic.
PHYS = {"Cstar": (100.0, 2500.0), "Cset": (300.0, 2500.0)}


def physical_ok(names, beta) -> bool:
    for n, v in zip(names, beta):
        if n in PHYS and not (PHYS[n][0] <= v <= PHYS[n][1]):
            return False
    return True


# ------------------------------------------------------------------ data
def load_system(name: str, end: str | None = None, start: str | None = None):
    w = pd.read_csv(CLEAN / "water_quality_daily.csv", parse_dates=["date"])
    g = w[w.system == name].sort_values("date").reset_index(drop=True)
    if end is not None:
        g = g[g.date <= pd.Timestamp(end)].reset_index(drop=True)
    if start is not None:
        g = g[g.date >= pd.Timestamp(start)].reset_index(drop=True)
    # readings flagged as isolated spikes in the audit are left out of the fit (not deleted)
    fl = pd.read_csv(CLEAN / "flagged_readings.csv", parse_dates=["date"])
    fl = fl[(fl.system == name) & (fl.parameter == "conductivity")]
    g = g[~g.date.isin(fl.date)].reset_index(drop=True)
    days = pd.date_range(g.date.min(), g.date.max(), freq="D")
    temp = g.set_index("date").temperature.reindex(days).interpolate().bfill().ffill().values
    return g, days, temp


def dosing_windows(g: pd.DataFrame):
    """Periods where conductivity is held almost constant -> automatic dosing assumed ON.
    Inferred from the data (facility manager: both Main and Tecniplast have dosing that works only
    part of the time)."""
    s = g.set_index("date").conductivity
    sd = s.rolling(STABLE_MIN_READINGS, center=True, min_periods=4).std()
    on = sd < STABLE_SD
    runs = (on != on.shift()).cumsum()
    return [(r.index[0], r.index[-1]) for _, r in on.groupby(runs)
            if r.iloc[0] and len(r) >= STABLE_MIN_READINGS]


# ------------------------------------------------------------------ ODE building blocks
def propagate(lam, forcing, impulse):
    """Exact daily solution of dC/dt = -lam_i*C + forcing_i on day i, with an instantaneous
    jump impulse_i at the start of day i. Returns C at the start of each day (after jump)."""
    n = len(lam)
    out = np.empty(n)
    c = 0.0
    e = np.exp(-lam)
    gain = (1 - e) / lam
    for i in range(n):
        c += impulse[i]
        out[i] = c
        c = c * e[i] + forcing[i] * gain[i]
    return out


class ConductivityModel:
    def __init__(self, name: str, end: str | None = None, start: str | None = None):
        self.name = name
        self.g, self.days, self.temp = load_system(name, end, start)
        self.nd = len(self.days)
        self.obs = ((self.g.date - self.days[0]).dt.days).values
        self.y = self.g.conductivity.values.astype(float)
        self.windows = dosing_windows(self.g)
        self.active = np.zeros(self.nd)
        for a, b in self.windows:
            self.active[(self.days >= a) & (self.days <= b)] = 1.0
        self.has_ctrl = bool(self.active.any())
        self.tdev = self.temp - self.temp.mean()
        # candidate jump days: middle day of each gap between readings, outside dosing windows
        cand = []
        for a, b in zip(self.obs[:-1], self.obs[1:]):
            if self.active[a] and self.active[b]:
                continue
            cand.append(a + max(1, (b - a) // 2))
        self.cand = np.array(sorted(set(cand)))

    def blocks(self, r, k, use_temp):
        lam = r + k * self.active
        z = np.zeros(self.nd)
        imp0 = z.copy(); imp0[0] = 1.0
        cols = {"C0": propagate(lam, z, imp0),
                "Cstar": propagate(lam, r * np.ones(self.nd), z)}
        if use_temp:
            cols["g"] = propagate(lam, r * self.tdev, z)
        if self.has_ctrl:
            cols["Cset"] = propagate(lam, k * self.active, z)
        base = np.column_stack([v[self.obs] for v in cols.values()])
        jumps = np.empty((len(self.obs), len(self.cand)))
        for j, d in enumerate(self.cand):
            imp = z.copy(); imp[d] = 1.0
            jumps[:, j] = propagate(lam, z, imp)[self.obs]
        return list(cols), base, jumps

    def fit_given(self, r, k, use_temp):
        """Best jump set for fixed r, k: LASSO path for selection, OLS refit, BIC choice."""
        names, B, J = self.blocks(r, k, use_temp)
        n = len(self.y)
        # remove the unpenalised part (Frisch-Waugh) so only jumps are penalised
        Q, _ = np.linalg.qr(B)
        P = lambda M: M - Q @ (Q.T @ M)
        yt, Jt = P(self.y), P(J)
        scale = np.linalg.norm(Jt, axis=0); scale[scale == 0] = 1
        _, coefs, _ = lasso_path(Jt / scale, yt, n_alphas=60, eps=1e-3)
        best = None
        seen = set()
        for c in coefs.T:
            sup = tuple(np.flatnonzero(np.abs(c) > 1e-10))
            if sup in seen or len(sup) > n // 3:
                continue
            seen.add(sup)
            X = np.column_stack([B, J[:, list(sup)]]) if sup else B
            beta, *_ = np.linalg.lstsq(X, self.y, rcond=None)
            D = beta[B.shape[1]:]
            keep = [s for s, dj in zip(sup, D) if abs(dj) >= MIN_JUMP]
            if len(keep) != len(sup):  # drop tiny jumps and refit
                X = np.column_stack([B, J[:, keep]]) if keep else B
                beta, *_ = np.linalg.lstsq(X, self.y, rcond=None)
                sup = tuple(keep)
            if not physical_ok(names, beta):
                continue
            rss = float(np.sum((X @ beta - self.y) ** 2))
            npar = B.shape[1] + 1 + (1 if self.has_ctrl else 0)  # + r (+ k)
            bic = n * np.log(rss / n) + (npar + JUMP_PENALTY * len(sup)) * np.log(n)
            if best is None or bic < best["bic"]:
                best = {"bic": bic, "rss": rss, "beta": beta, "sup": sup, "names": names,
                        "r": r, "k": k if self.has_ctrl else 0.0, "use_temp": use_temp}
        return best

    def fit(self, r_grid=R_GRID, use_temp=False):
        best = None
        for r in r_grid:
            for k in (K_GRID if self.has_ctrl else [0.0]):
                f = self.fit_given(float(r), float(k), use_temp)
                if f is not None and (best is None or f["bic"] < best["bic"]):
                    best = f
        self.best = best
        nb = len(best["names"])
        self.params = dict(zip(best["names"], best["beta"][:nb]))
        self.params.update({"r": best["r"], "k": best["k"]})
        self.jump_days = self.cand[list(best["sup"])]
        self.jumps = best["beta"][nb:]
        return self

    # ---------------------------------------------------------- simulation
    def simulate_daily(self):
        p = self.params
        lam = p["r"] + p["k"] * self.active
        forcing = p["r"] * (p["Cstar"] + p.get("g", 0.0) * self.tdev) + p["k"] * self.active * p.get("Cset", 0.0)
        imp = np.zeros(self.nd); imp[0] = p["C0"]
        for d, D in zip(self.jump_days, self.jumps):
            imp[d] += D
        return propagate(lam, forcing, imp)

    def simulate_ivp(self, per_day=24):
        p = self.params
        c = 0.0
        jumps = np.zeros(self.nd); jumps[0] = p["C0"]
        for d, D in zip(self.jump_days, self.jumps):
            jumps[d] += D
        tt, cc = [], []
        for i in range(self.nd):
            c += jumps[i]
            target = p["Cstar"] + p.get("g", 0.0) * self.tdev[i]
            a = self.active[i]
            rhs = lambda t, x: p["r"] * (target - x) + a * p["k"] * (p.get("Cset", 0.0) - x)
            s = solve_ivp(rhs, (i, i + 1), [c], t_eval=np.linspace(i, i + 1, per_day + 1), rtol=1e-9, atol=1e-9)
            tt.append(s.t[:-1]); cc.append(s.y[0][:-1]); c = s.y[0][-1]
        return np.concatenate(tt), np.concatenate(cc)


def metrics(y, yhat):
    e = yhat - y
    return {"RMSE": float(np.sqrt(np.mean(e ** 2))), "MAE": float(np.mean(np.abs(e))),
            "R2": float(1 - np.sum(e ** 2) / np.sum((y - y.mean()) ** 2))}


# ------------------------------------------------------------------ main
def main():
    vol = pd.read_csv(CLEAN / "systems.csv").set_index("system").volume_l
    out = ["# Step 2 - Conductivity mass-balance model\n",
           "Fitted on ALL readings of each system, to show the model structure and the inferred",
           "events. The calibration (Apr-Aug) / validation (Sep-Oct) split and the parameter",
           "uncertainty come in Step 4.\n"]
    prow, erow = [], []
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.8), sharex=True)
    for ax, name in zip(axes, SYSTEMS):
        m = ConductivityModel(name).fit()
        bic0 = m.best["bic"]
        mt = ConductivityModel(name).fit(use_temp=True)
        mr = ConductivityModel(name).fit(r_grid=[0.10])
        d_temp, d_r10 = mt.best["bic"] - bic0, mr.best["bic"] - bic0
        # accept the temperature term only if clearly better AND identifiable: when r sits at the
        # edge of its range, r and the temperature effect cannot be told apart (g blows up)
        r_edge = mt.params["r"] in (R_GRID[0], R_GRID[-1])
        temp_ok = d_temp < -2 and not r_edge
        if temp_ok:
            m = mt
        p = m.params
        pred = m.simulate_daily()
        t_ivp, c_ivp = m.simulate_ivp()
        check = float(np.max(np.abs(c_ivp[::24] - pred)))
        met = metrics(m.y, pred[m.obs])
        up, down = int((m.jumps > 0).sum()), int((m.jumps < 0).sum())

        out += [f"## {name}  (volume {vol[name]:.0f} L)\n",
                f"- Readings used: {len(m.y)} ({m.days[0].date()} to {m.days[-1].date()})",
                "- Dosing controller ON (inferred from flat stretches): " +
                (", ".join(f"{a.date()} to {b.date()}" for a, b in m.windows) or "none detected"),
                f"- Effective renewal rate r = {p['r']:.3f} per day (facility assumption 0.10); "
                f"a jump loses half its size in {np.log(2)/p['r']:.1f} days",
                f"- Forcing r = 0.10 instead: BIC {d_r10:+.1f} "
                + ("(clearly worse)" if d_r10 > 10 else "(worse)" if d_r10 > 2 else "(about as good)"),
                f"- Equilibrium without intervention C* = {p['Cstar']:.0f} µS/cm"]
        if m.has_ctrl:
            out.append(f"- Controller set-point C_set = {p['Cset']:.0f} µS/cm (strength k = {p['k']:.1f} per day)")
        out += [f"- Temperature effect on evaporation: BIC {d_temp:+.1f} -> "
                + ("supported, included" if temp_ok else
                   "better BIC but not identifiable (renewal rate at the edge of its range), left out" if d_temp < -2
                   else "not supported by the data, left out")
                + (f" (g = {p['g']:.1f} µS/cm per °C)" if "g" in p else ""),
                f"- Inferred jumps: {len(m.jumps)} ({up} up = salt additions, {down} down = dilutions)",
                f"- In-sample fit: RMSE {met['RMSE']:.1f} µS/cm, MAE {met['MAE']:.1f}, R² {met['R2']:.3f}",
                f"- Check: exact daily solution vs solve_ivp differ by at most {check:.1e} µS/cm\n"]
        prow.append({"system": name, **p, "n_jumps": len(m.jumps), "dBIC_temperature": d_temp,
                     "dBIC_r_fixed_0.10": d_r10, **met})
        for dd, D in zip(m.jump_days, m.jumps):
            erow.append({"system": name, "date": m.days[dd].date(), "jump_uS_cm": round(float(D), 1),
                         "type": "salt addition (inferred)" if D > 0 else "dilution / water change (inferred)",
                         "approx_salt_g": round(float(D) * vol[name] * MG_PER_L_PER_US / 1000, 1) if D > 0 else None})

        for a, b in m.windows:
            ax.axvspan(a, b, color="#9a9994", alpha=0.18, lw=0)
        ax.plot(m.days[0] + pd.to_timedelta(t_ivp, unit="D"), c_ivp, color="#52514e", lw=1.2, label="model")
        ax.plot(m.g.date, m.y, "o", ms=3.5, color=COLORS[name], label="measured")
        for dd, D in zip(m.jump_days, m.jumps):
            y0 = pred[dd]
            ax.annotate("", xy=(m.days[dd], y0), xytext=(m.days[dd], y0 - D),
                        arrowprops=dict(arrowstyle="-|>", color="#0b0b0b", lw=0.8, shrinkA=0, shrinkB=0))
        ax.set_title(f"{name}  ·  r = {p['r']:.2f}/day  ·  {len(m.jumps)} inferred jumps  ·  "
                     f"RMSE {met['RMSE']:.0f} µS/cm", loc="left", fontsize=9, color="#52514e")
        ax.set_ylabel("µS/cm")
        if name == SYSTEMS[0]:
            ax.legend(frameon=False, loc="upper right", fontsize=8)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.suptitle("Conductivity: mass-balance model vs measurements\n(arrows = inferred jumps; "
                 "grey bands = automatic dosing ON, inferred)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "step2_conductivity_model.png", dpi=150)
    plt.close(fig)

    pd.DataFrame(prow).to_csv(CLEAN / "step2_conductivity_params.csv", index=False)
    ev = pd.DataFrame(erow)
    ev.to_csv(CLEAN / "step2_inferred_jumps.csv", index=False)
    out += ["## Inferred jumps (NOT observed: estimated from the data)\n", ev.to_string(index=False),
            f"\napprox_salt_g: rough conversion {MG_PER_L_PER_US} mg/L per µS/cm (NaCl-like salt) x system "
            "volume. An order of magnitude, not a measurement."]
    (REP / "step2_conductivity.md").write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
