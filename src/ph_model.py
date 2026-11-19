"""Step 3 - Process-based (mechanistic) model of pH, one per water system.

PLAIN-LANGUAGE IDEA
Fish breathe out CO2 and their waste is turned into acid by the biofilter bacteria
(nitrification). Both push pH DOWN, and more fish per litre push harder. Renewal
replaces part of the water every day with buffered makeup water, which pulls pH back
towards the makeup pH. Automatic dosing (when it works) and cleanings can also move pH.

EQUATION (ODE; pH in pH units, t in days)

    dpH/dt = r * (pH_in - pH)                   renewal pulls pH towards the makeup pH
           - alpha * F(t) / V                   acid/CO2 load: F = fish in the system, V = volume (L)
           + a(t) * k * (pH_set - pH)           automatic dosing ON (a = 1) in the dosing windows
           + delta * c(t)                       cleaning: c(t) = fraction of tanks cleaned that day
           + sum_j D_j * delta(t - t_j)         other jumps (e.g. unlogged buffer additions), inferred

Parameters, in plain words:
  r       effective renewal rate per day (how fast pH returns to the makeup value)
  pH_in   pH the system would reach with no fish (makeup pH, as seen by the system)
  alpha   how much pH each fish-per-litre removes per day; at equilibrium pH = pH_in - alpha*F/(r*V)
  pH_set  dosing set-point; k = how strongly the controller acts
  delta   pH change caused by cleaning all tanks of a system on one day

SIMPLIFICATION (stated openly): pH is a logarithmic quantity and real carbonate chemistry
is non-linear. Modelling pH itself with linear terms is a first-order approximation that
holds for small changes around the operating point (pH 7-8.5) and constant alkalinity.
We have no alkalinity or CO2 data, so a full carbonate model cannot be parameterised.

FITTING: identical machinery to Step 2 (exact daily ODE solution, grid over r and k, LASSO
for which gaps need a jump, OLS refit, modified BIC). Each optional term (fish load,
cleaning, dosing) is kept only if it lowers BIC by more than 2. The final model is checked
against scipy.integrate.solve_ivp.

Run: python src/ph_model.py
"""
from __future__ import annotations

import itertools
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

from conductivity_model import dosing_windows, propagate

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CLEAN, FIG, REP = ROOT / "data" / "clean", ROOT / "figures", ROOT / "reports"
SYSTEMS = ["Main System", "SA Tecniplast", "SA Aquaneering"]
COLORS = {"Main System": "#2a78d6", "SA Tecniplast": "#eb6834", "SA Aquaneering": "#1baf7a"}
TANKS = {"Main System": 53, "SA Tecniplast": 11, "SA Aquaneering": 70}  # typical tank count (audit)
MIN_JUMP = 0.25        # pH units; smaller changes treated as noise (meter resolution is 0.1)
JUMP_PENALTY = 3.0     # modified BIC, as in Step 2
R_GRID = np.round(np.geomspace(0.01, 1.0, 30), 4)
K_GRID = [0.5, 2.0, 10.0]
# Physically plausible ranges (ASSUMPTION): makeup/equilibrium pH and set-point between 5.5 and 9.5
PHYS = {"pH_in": (5.5, 9.5), "pH_set": (5.5, 9.5)}


def physical_ok(names, beta) -> bool:
    for n, v in zip(names, beta):
        if n in PHYS and not (PHYS[n][0] <= v <= PHYS[n][1]):
            return False
    return True


def load(name: str, end: str | None = None, start_override: str | None = None):
    sysinfo = pd.read_csv(CLEAN / "systems.csv").set_index("system")
    w = pd.read_csv(CLEAN / "water_quality_daily.csv", parse_dates=["date"])
    g = w[w.system == name].sort_values("date")
    start = pd.Timestamp(sysinfo.loc[name, "fish_reliable_from"])  # fish count must be known
    if start_override is not None:
        start = max(start, pd.Timestamp(start_override))
    g = g[g.date >= start]
    if end is not None:
        g = g[g.date <= pd.Timestamp(end)]
    fl = pd.read_csv(CLEAN / "flagged_readings.csv", parse_dates=["date"])
    fl = fl[(fl.system == name) & (fl.parameter == "ph")]
    g_fit = g[~g.date.isin(fl.date)].reset_index(drop=True)
    days = pd.date_range(g_fit.date.min(), g_fit.date.max(), freq="D")
    fish = pd.read_csv(CLEAN / "fish_per_system_daily.csv", parse_dates=["date"]).set_index("date")[name]
    fish = fish.reindex(days).ffill().bfill().values
    cl = pd.read_csv(CLEAN / "cleanings_per_system_daily.csv", parse_dates=["date"])
    cl = cl[cl.system == name].set_index("date").cleanings.reindex(days).fillna(0).values
    clean_frac = np.minimum(1.0, cl / TANKS[name])
    # dosing windows come from the CONDUCTIVITY record (flat stretches = controller working);
    # whether the same controller also holds pH is tested, not assumed
    cw = w[w.system == name].sort_values("date")
    windows = dosing_windows(cw)
    vol = float(sysinfo.loc[name, "volume_l"])
    return g_fit, days, fish / vol, clean_frac, windows, vol


class PHModel:
    def __init__(self, name, end=None, start=None):
        self.name = name
        self.g, self.days, self.fpl, self.clean, self.windows, self.vol = load(name, end, start)
        self.nd = len(self.days)
        self.obs = ((self.g.date - self.days[0]).dt.days).values
        self.y = self.g.ph.values.astype(float)
        self.active = np.zeros(self.nd)
        for a, b in self.windows:
            self.active[(self.days >= a) & (self.days <= b)] = 1.0
        # pH relaxes within days, so a jump's estimated size depends strongly on where in a gap it
        # is placed. We place it on the day of the next reading: the SMALLEST jump consistent with
        # the data (conservative choice, stated as an assumption).
        cand = [b for a, b in zip(self.obs[:-1], self.obs[1:])]
        self.cand = np.array(sorted(set(cand)))

    def blocks(self, r, k, terms):
        lam = r + k * self.active * ("dosing" in terms)
        z = np.zeros(self.nd)
        imp0 = z.copy(); imp0[0] = 1
        cols = {"pH0": propagate(lam, z, imp0), "pH_in": propagate(lam, r * np.ones(self.nd), z)}
        if "fish" in terms:
            cols["alpha"] = propagate(lam, -self.fpl, z)
        if "cleaning" in terms:
            cols["delta"] = propagate(lam, z, self.clean)
        if "dosing" in terms:
            cols["pH_set"] = propagate(lam, k * self.active, z)
        B = np.column_stack([v[self.obs] for v in cols.values()])
        J = np.empty((len(self.obs), len(self.cand)))
        for j, d in enumerate(self.cand):
            imp = z.copy(); imp[d] = 1
            J[:, j] = propagate(lam, z, imp)[self.obs]
        return list(cols), B, J

    def fit_given(self, r, k, terms):
        names, B, J = self.blocks(r, k, terms)
        n = len(self.y)
        Q, _ = np.linalg.qr(B)
        yt, Jt = self.y - Q @ (Q.T @ self.y), J - Q @ (Q.T @ J)
        sc = np.linalg.norm(Jt, axis=0); sc[sc == 0] = 1
        _, coefs, _ = lasso_path(Jt / sc, yt, n_alphas=50, eps=1e-3)
        best, seen = None, set()
        for c in np.column_stack([np.zeros(len(self.cand)), coefs]).T:
            sup = tuple(np.flatnonzero(np.abs(c) > 1e-12))
            if sup in seen or len(sup) > n // 4:
                continue
            seen.add(sup)
            X = np.column_stack([B, J[:, list(sup)]]) if sup else B
            beta, *_ = np.linalg.lstsq(X, self.y, rcond=None)
            keep = [s for s, dj in zip(sup, beta[B.shape[1]:]) if abs(dj) >= MIN_JUMP]
            if len(keep) != len(sup):
                X = np.column_stack([B, J[:, keep]]) if keep else B
                beta, *_ = np.linalg.lstsq(X, self.y, rcond=None)
                sup = tuple(keep)
            if not physical_ok(names, beta):
                continue
            rss = float(np.sum((X @ beta - self.y) ** 2))
            npar = B.shape[1] + 1 + ("dosing" in terms)
            bic = n * np.log(rss / n) + (npar + JUMP_PENALTY * len(sup)) * np.log(n)
            if best is None or bic < best["bic"]:
                best = dict(bic=bic, beta=beta, sup=sup, names=names, r=r,
                            k=k if "dosing" in terms else 0.0, terms=terms)
        return best

    def fit(self, terms, r_grid=R_GRID):
        best = None
        for r in r_grid:
            for k in (K_GRID if "dosing" in terms else [0.0]):
                f = self.fit_given(float(r), float(k), terms)
                if f is not None and (best is None or f["bic"] < best["bic"]):
                    best = f
        return best

    def set_fit(self, best):
        self.best = best
        nb = len(best["names"])
        self.p = dict(zip(best["names"], best["beta"][:nb]))
        self.p.update(r=best["r"], k=best["k"])
        self.jump_days = self.cand[list(best["sup"])]
        self.jumps = best["beta"][nb:]
        return self

    def simulate(self):
        p = self.p
        act = self.active * ("pH_set" in p)
        lam = p["r"] + p["k"] * act
        forcing = p["r"] * p["pH_in"] - p.get("alpha", 0) * self.fpl + p["k"] * act * p.get("pH_set", 0)
        imp = p.get("delta", 0) * self.clean
        imp[0] += p["pH0"]
        for d, D in zip(self.jump_days, self.jumps):
            imp[d] += D
        return propagate(lam, forcing, imp)

    def simulate_ivp(self, per_day=24):
        p = self.p
        act = self.active * ("pH_set" in p)
        imp = p.get("delta", 0) * self.clean
        imp[0] += p["pH0"]
        for d, D in zip(self.jump_days, self.jumps):
            imp[d] += D
        x, tt, cc = 0.0, [], []
        for i in range(self.nd):
            x += imp[i]
            rhs = lambda t, v, i=i: (p["r"] * (p["pH_in"] - v) - p.get("alpha", 0) * self.fpl[i]
                                     + act[i] * p["k"] * (p.get("pH_set", 0) - v))
            s = solve_ivp(rhs, (i, i + 1), [x], t_eval=np.linspace(i, i + 1, per_day + 1), rtol=1e-10, atol=1e-10)
            tt.append(s.t[:-1]); cc.append(s.y[0][:-1]); x = s.y[0][-1]
        return np.concatenate(tt), np.concatenate(cc)


def select_terms(m: PHModel):
    """Fit every combination of optional terms; keep the lowest BIC. Report, for each term,
    how much BIC worsens if it is removed from the best model (positive = term is useful)."""
    options = ["fish", "cleaning"] + (["dosing"] if m.active.any() else [])
    fits = {}
    for nterms in range(len(options) + 1):
        for combo in itertools.combinations(options, nterms):
            f = m.fit(set(combo))
            fits[frozenset(combo)] = f if f is not None else {"bic": np.inf}
    best_set = min(fits, key=lambda s: fits[s]["bic"])
    # prefer the simpler model unless the extra term improves BIC by > 2
    for s in sorted(fits, key=len):
        if fits[s]["bic"] < fits[best_set]["bic"] + 2 and s <= best_set:
            best_set = s
            break
    contrib = {}  # BIC(with term) - BIC(without term), around the chosen model; negative = term helps
    for t in options:
        with_t, without_t = best_set | {t}, best_set - {t}
        contrib[t] = fits[frozenset(with_t)]["bic"] - fits[frozenset(without_t)]["bic"]
    return fits[best_set], best_set, contrib, options


def metrics(y, yhat):
    e = yhat - y
    return {"RMSE": float(np.sqrt(np.mean(e ** 2))), "MAE": float(np.mean(np.abs(e))),
            "R2": float(1 - np.sum(e ** 2) / np.sum((y - y.mean()) ** 2))}


def main():
    out = ["# Step 3 - pH model\n",
           "Fitted on all readings from the date each system's fish count is reliable. The",
           "calibration/validation split and uncertainty come in Step 4.\n"]
    rows, erows = [], []
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.8), sharex=True)
    for ax, name in zip(axes, SYSTEMS):
        m = PHModel(name)
        best, terms, contrib, options = select_terms(m)
        m.set_fit(best)
        p = m.p
        pred = m.simulate()
        t_ivp, c_ivp = m.simulate_ivp()
        check = float(np.max(np.abs(c_ivp[::24] - pred)))
        met = metrics(m.y, pred[m.obs])
        fpl = m.fpl[m.obs]
        out += [f"## {name}  (volume {m.vol:.0f} L, {fpl.min():.1f}-{fpl.max():.1f} fish per litre)\n",
                f"- Readings used: {len(m.y)} ({m.days[0].date()} to {m.days[-1].date()})",
                f"- Terms kept: {', '.join(sorted(terms)) or 'renewal only'}",
                "- BIC with vs without each term (negative = the term helps; kept only if < -2): "
                + ", ".join(f"{t} {contrib[t]:+.1f}" for t in options),
                f"- pH relaxation rate r = {p['r']:.3f} per day (half-return time {np.log(2)/p['r']:.1f} days). "
                "For pH this includes CO2 degassing by aeration, not only water renewal."]
        if "alpha" in p:
            sens = p["alpha"] / p["r"]
            eq = p["pH_in"] - sens * fpl.mean()
            out += [f"- Fish load alpha = {p['alpha']:.4f} pH per day per fish/L; equilibrium sensitivity "
                    f"alpha/r = {sens:.2f} pH units lower per extra fish per litre",
                    f"- Equilibrium pH at the mean density ({fpl.mean():.1f} fish/L) = {eq:.2f}",
                    f"- pH_in (= no-fish value) = {p['pH_in']:.2f} is an EXTRAPOLATION far outside the observed "
                    f"density range ({fpl.min():.1f}-{fpl.max():.1f} fish/L) and is not physically meaningful here",
                    "- CAUTION: fish numbers only went down steadily in this period, so the fish term cannot be "
                    "separated from any other slow drift (biofilter maturation, season, routine changes)."
                    if np.all(np.diff(fpl) <= 1e-9) else "- Fish numbers changed in steps, which helps identify the effect."]
        else:
            out.append(f"- Equilibrium pH (no detectable fish effect) = {p['pH_in']:.2f}")
        if "delta" in p:
            out.append(f"- Cleaning effect delta = {p['delta']:+.2f} pH if all tanks are cleaned on one day")
        if "pH_set" in p:
            out.append(f"- Dosing set-point pH_set = {p['pH_set']:.2f} (k = {p['k']:.1f} per day)")
        out += [f"- Other inferred jumps: {len(m.jumps)}",
                f"- In-sample fit: RMSE {met['RMSE']:.3f}, MAE {met['MAE']:.3f}, R² {met['R2']:.3f}",
                f"- Check vs solve_ivp: max difference {check:.1e}\n"]
        rows.append({"system": name, **p, "terms": "+".join(sorted(terms)), "n_jumps": len(m.jumps),
                     **{f"dBIC_{t}": contrib[t] for t in options}, **met})
        for d, D in zip(m.jump_days, m.jumps):
            erows.append({"system": name, "date": m.days[d].date(), "jump_pH": round(float(D), 2)})

        for a, b in m.windows:
            ax.axvspan(a, b, color="#9a9994", alpha=0.18, lw=0)
        ax.plot(m.days[0] + pd.to_timedelta(t_ivp, unit="D"), c_ivp, color="#52514e", lw=1.2, label="model")
        ax.plot(m.g.date, m.y, "o", ms=3.5, color=COLORS[name], label="measured")
        for i in np.flatnonzero(m.clean > 0.2):
            ax.axvline(m.days[i], color="#0b0b0b", lw=0.8, ls=":")
        for d, D in zip(m.jump_days, m.jumps):
            ax.annotate("", xy=(m.days[d], pred[d]), xytext=(m.days[d], pred[d] - D),
                        arrowprops=dict(arrowstyle="-|>", color="#0b0b0b", lw=0.8, shrinkA=0, shrinkB=0))
        ax.set_title(f"{name}  ·  terms: {', '.join(sorted(terms)) or 'renewal only'}  ·  "
                     f"RMSE {met['RMSE']:.2f}  ·  R² {met['R2']:.2f}", loc="left", fontsize=9, color="#52514e")
        ax.set_ylabel("pH")
        if name == SYSTEMS[0]:
            ax.legend(frameon=False, loc="lower left", fontsize=8)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.suptitle("pH: process model vs measurements\n(dotted = big cleaning days; arrows = other inferred "
                 "jumps; grey = dosing ON per conductivity record)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "step3_ph_model.png", dpi=150)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(CLEAN / "step3_ph_params.csv", index=False)
    ev = pd.DataFrame(erows)
    ev.to_csv(CLEAN / "step3_ph_jumps.csv", index=False)
    out += ["## Other inferred pH jumps (NOT observed)\n", ev.to_string(index=False) if len(ev) else "None."]
    (REP / "step3_ph.md").write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
