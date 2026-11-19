---
title: "Process-based and machine-learning models of conductivity and pH in three zebrafish recirculating systems: an exploratory study with routine facility records"
author: "Bruno Ferreira"
date: "November 2026"
---

# Abstract

Routine water-quality logs from research aquaculture facilities are rarely used for modelling. We used six months of manual daily readings (temperature, conductivity, pH) from three independent zebrafish (*Danio rerio*) recirculating systems (Apr-Nov 2026; 74-134 readings per system), together with tank-level records of fish numbers and cleanings, to build and test (i) a salt mass-balance ordinary differential equation (ODE) for conductivity, (ii) a linearised pH model with fish load, dosing and cleaning terms, and (iii) small machine-learning (ML) forecasters. Unlogged salt additions and dilutions were inferred as sparse jumps using LASSO selection and a modified BIC; intermittent automatic dosing was identified from flat stretches of the record. In-sample, the conductivity model explained 81-95 % of the variance (RMSE 37-55 µS/cm) and showed effective renewal rates of 0.08, 0.30 and 0.03 d⁻¹, all statistically inconsistent with the nominal 10 %/day; the implied absolute exchange was similar for two systems (≈ 42-52 L/day) and much lower for the third (≈ 4 L/day). With a strict time split (calibration Apr-Aug, validation Sep-Oct), neither the mechanistic models nor ridge regression or gradient boosting beat a persistence forecast for conductivity, because the dominant drivers were unlogged human interventions. For pH in one system, which relaxes to its set level within about two days, ML models and the ODE beat persistence (best RMSE 0.318 vs 0.370 pH units). We discuss what additional logging would make such models useful for predictive monitoring.

# 1. Introduction

Predictive monitoring of water quality is a building block of early-warning systems for disease in aquaculture, where parasite outbreaks are linked to environmental conditions. Two modelling traditions compete: process-based (mechanistic) models, which encode mass balances and have interpretable parameters, and data-driven ML models, which need fewer assumptions but more data. Comparing them honestly on small, imperfect operational datasets is an important but rarely reported exercise.

This exploratory personal project uses records already collected in a zebrafish facility through a facility-management app. The aims were: (1) audit the data; (2) build process-based models of conductivity and pH; (3) estimate their parameters with uncertainty, using a calibration/validation time split; (4) test whether parameters transfer to an independent system; and (5) compare against simple ML forecasters under the same split.

# 2. Data

**Systems.** Three independent recirculating systems: *Main System* (3 × 175 L = 525 L), *SA Tecniplast* (175 L) and *SA Aquaneering* (150 L). Main and Tecniplast were logged from 23 Apr 2026, Aquaneering from 8 Jul 2026 (earlier readings were excluded at the facility manager's request); data end on 10 Nov 2026. Makeup water is reverse-osmosis water with small additions of salt and/or pH buffer. Main and Tecniplast have automatic dosing that worked only part of the time.

**Water quality.** Manual readings of temperature (°C), conductivity (µS/cm) and pH, entered in the app (369 rows; 4 deleted entries and 2 duplicated days removed; 28 readings from two short-lived systems excluded). Coverage was 59-66 % of calendar days (weekends and a mid-August to mid-September gap). Although the protocol specifies 09:00, the median recording time was ≈ 11:00 (only 4 % between 08:30 and 09:30); readings are treated as daily late-morning snapshots. Units were inferred from value ranges (no unit field). Four isolated spikes (a reading far from both neighbours, which agree with each other) were flagged and excluded from fitting but not deleted.

**Tank records.** 243 tank records with setup date, archive date, rack (mapped to systems), fish count and a dated history of moves and mortality. Daily fish per system was reconstructed backwards from mortality events. For Main, counts before 22 May 2026 reflect tanks still being entered into the app and were set to unknown. Fish density was 2.4-2.6 (Main), 0.7-1.9 (Tecniplast; fish added on 14 Sep) and 5.1-7.1 fish/L (Aquaneering, declining through mortality). Tank cleanings were logged mainly from August onwards. Feeding logs exist but were judged unreliable by the facility manager and were not used; water changes and manual salt additions are not logged.

**Anonymisation.** Published data contain no personal information (user e-mails, names, notes and project/genetic-line labels removed).

# 3. Methods

## 3.1 Conductivity: salt mass balance

$$\frac{dC}{dt} = r\,\big(C^{*} + g\,(T-\bar T) - C\big) \;+\; a(t)\,k\,(C_{set}-C) \;+\; \sum_j D_j\,\delta(t-t_j)$$

*C* is conductivity; *r* the effective renewal rate (fraction of volume exchanged per day); *C\** the equilibrium conductivity reached without intervention, which combines makeup-water conductivity $C_{in}$ and evaporative concentration, $C^{*}=C_{in}(1+E/rV)$ (*E*, evaporation in L/day; *V*, volume), $C_{in}$ and *E* cannot be separated without measuring the makeup water; *g* an optional temperature effect on evaporation; *a(t)* = 1 during automatic-dosing periods, with set-point $C_{set}$ and strength *k*; $D_j$ are jumps at times $t_j$: inferred manual salt additions (*D* > 0) or dilutions/water changes (*D* < 0). Dosing periods were identified as runs of ≥ 5 readings with rolling SD < 15 µS/cm.

## 3.2 pH: linearised acid/CO₂ balance

$$\frac{dpH}{dt} = r\,(pH_{in}-pH) \;-\; \alpha\,\frac{F(t)}{V} \;+\; a(t)\,k\,(pH_{set}-pH) \;+\; \delta_c\,c(t) \;+\; \sum_j D_j\,\delta(t-t_j)$$

*r* is a relaxation rate (renewal plus CO₂ degassing by aeration); $pH_{in}$ the no-fish equilibrium; *α* the acidifying effect of fish density *F/V* (CO₂ and nitrification); *c(t)* the fraction of tanks cleaned on a day and $\delta_c$ its effect; dosing and jumps as above. Modelling pH directly with linear terms is a first-order approximation of carbonate chemistry, valid for small deviations at constant alkalinity; no alkalinity or CO₂ data were available.

## 3.3 Estimation

Both ODEs are linear in *C* (or pH) with piecewise-constant daily inputs, so they were solved exactly day by day; the final fits were re-simulated with `scipy.integrate.solve_ivp` (maximum difference < 10⁻⁶). For fixed *r* and *k* (grid search), the model output is linear in the remaining parameters and in a candidate jump in every gap between readings. Jumps were selected with the LASSO path and refitted by ordinary least squares; the model with the lowest modified BIC (each jump penalised 3 ln *n*: size plus location) was retained, with a minimum jump of 75 µS/cm or 0.25 pH units. Optional terms (temperature, fish load, cleaning, dosing) were kept only if ΔBIC < −2. Fits outside physically plausible ranges (*C\** 100-2500 µS/cm; $pH_{in}$, $pH_{set}$ 5.5-9.5) were rejected. Conductivity jumps are placed mid-gap; pH jumps on the day of the next reading (the smallest jump consistent with the data).

## 3.4 Calibration, validation and uncertainty

Models were calibrated on readings up to 31 Aug 2026 and validated on 1 Sep-31 Oct 2026 with next-reading forecasts: each forecast starts from the previous measured value and runs the calibrated ODE, with logged/measured inputs (temperature, fish, cleanings, dosing periods) treated as known and no knowledge of future jumps. Parameter uncertainty was estimated with a moving-block residual bootstrap (300 resamples, block = 5 readings) conditional on the inferred jump dates. Parameter stability was checked by refitting on Sep-Oct only.

## 3.5 Transfer test and ML comparison

For Aquaneering we compared (A) full transfer of donor parameters (*r* and equilibrium) from Main or Tecniplast, (B) transfer of *r* only, and (C) an own refit; jumps were always re-inferred. ML forecasters predicted the change to the next reading from lags (1-3 readings), recent differences, gap length, previous temperature, fish density and logged cleanings: ridge regression (penalty by time-series cross-validation) and gradient boosting (depth 2, 150 trees), per system and pooled over systems. All methods were scored on identical test points against persistence (next = last reading), with 95 % bootstrap intervals for RMSE differences.

# 4. Results

## 4.1 Conductivity

| | Main | Tecniplast | Aquaneering |
|---|---|---|---|
| Renewal *r* (d⁻¹) | 0.080 | 0.300 | 0.027 |
| ΔBIC if *r* fixed at 0.10 | +11.7 | +29.8 | +56.0 |
| *C\** (µS/cm) | 709 | 651 | 1202 |
| $C_{set}$ (µS/cm) | 747 | 744 | to |
| Inferred jumps (up / down) | 7 / 4 | 6 / 0 | 3 / 3 |
| In-sample RMSE (µS/cm); R² | 37; 0.90 | 39; 0.81 | 55; 0.95 |

: Conductivity model, fitted on the full record.

The fitted renewal rates differ strongly between systems and none matches the nominal 10 %/day. Converted to volumes, Main and Tecniplast exchange a similar ≈ 42 and ≈ 52 L/day, whereas Aquaneering exchanges ≈ 4 L/day, consistent with its higher, slowly drifting conductivity. A temperature effect was supported only for Main (+62 µS/cm per °C), which may reflect evaporation or incomplete temperature compensation of the meter. Inferred salt additions correspond to 5-150 g of salt (rough conversion 0.5 mg/L per µS/cm).

## 4.2 pH

The data supported few terms (Table 2). No fish-load effect was detectable in Main or Tecniplast, even though Tecniplast nearly tripled its density on 14 Sep. Tecniplast supported dosing (set-point 7.68) and a cleaning effect (−0.49 pH for a full cleaning, from a single event). Aquaneering pH was not explained beyond its mean (R² ≈ 0).

| | Main | Tecniplast | Aquaneering |
|---|---|---|---|
| Terms supported | relaxation only | dosing, cleaning | relaxation only |
| *r* (d⁻¹) | 0.022 | 0.329 | 0.127 |
| Equilibrium pH | 8.02 | 7.43 (set-point 7.68) | 7.76 |
| In-sample RMSE; R² | 0.19; 0.73 | 0.19; 0.48 | 0.44; 0.01 |

: pH model, fitted from the date fish counts are reliable.

## 4.3 Calibration, validation and uncertainty

Calibration on Apr-Aug gave wide or boundary-limited intervals: e.g. Main *r* = 0.015 d⁻¹ (95 % CI 0.005-0.036), Tecniplast *r* = 0.056 (0.051-0.077) with *C\** at the lower plausible limit, Aquaneering *C\** at the upper limit. Tecniplast was held flat by its dosing controller for most of the calibration period, leaving little information on free dynamics; refitting on Sep-Oct gave *r* = 0.34, outside the calibration interval. In validation, next-reading forecasts of the ODEs matched persistence for Main (conductivity RMSE 41.7 vs 41.7; pH 0.243 vs 0.243) and were worse or indistinguishable for the other systems (Fig. 2).

## 4.4 Transfer to an independent system

Full transfer of donor parameters to Aquaneering failed for conductivity (RMSE 118 and 605 µS/cm vs 55 for the own fit), because its equilibrium level differs. Transferring only Main's renewal rate gave an almost equivalent fit (RMSE 59), whereas Tecniplast's fast rate did not (168). No variant explained Aquaneering pH.

## 4.5 ML vs mechanistic vs persistence

For conductivity, no method was clearly better than persistence (one marginal exception: per-system ridge for Tecniplast, 113 vs 115 µS/cm). Per-system gradient boosting overfitted (Main: 82.7 vs 41.7). For Tecniplast pH, pooled ridge (0.318), pooled gradient boosting (0.332), per-system ridge (0.327) and the ODE (0.365) beat persistence (0.370), with bootstrap intervals excluding zero.

# 5. Discussion

The process-based models summarise the systems in a few interpretable numbers, renewal rates, equilibria, dosing set-points and the dates and sizes of unlogged interventions, and produced a testable operational finding: renewal behaves like a similar absolute flow in two systems and is much lower in the third, rather than a fixed 10 % of volume. This is the main value of the mechanistic approach here.

As forecasters, however, neither mechanistic nor ML models improved on persistence for conductivity. The reason is structural: the largest changes are caused by unlogged human actions (manual salt additions, water changes, the dosing controller switching on and off), which no model can anticipate from past readings. Where the dynamics are autonomous and fast, pH in Tecniplast, which relaxes to its set level within about two days, both approaches beat persistence, and simple pooled linear models were the most reliable ML option. With ≈ 200 training rows in total, flexible ML overfitted unless data were pooled.

The parameter instability between periods is informative rather than merely a failure: it shows that the systems were operated differently over time (dosing working or not), which a time-split validation reveals and an in-sample fit hides.

**Practical recommendations** for turning routine logs into predictive tools: log every salt addition, buffer addition and water change (date, amount); measure makeup-water conductivity and pH periodically; record dosing-controller status; keep the reading time fixed; and add alkalinity measurements if pH is to be modelled mechanistically.

# 6. Limitations

- Exploratory, single-facility, six-month dataset; 74-134 readings per system with 34-41 % of days missing.
- All salt additions, dilutions and pH jumps are **inferred**, not observed, and their timing within gaps is assumed. Confidence intervals are conditional on the inferred jump dates.
- Dosing periods are inferred from flat stretches of the record, not from controller logs.
- Makeup composition is assumed constant; $C_{in}$ and evaporation cannot be separated.
- The pH model is a linear approximation without carbonate chemistry.
- Fish counts are reconstructed from app records and assume that tank moves and archiving were logged on time; feeding was not usable.
- Forecast tests have 32-46 points per system; differences of a few percent are not meaningful.

# 7. Reproducibility

All code is available in the accompanying repository. The facility data are not public; they belong to the host institution and are available from the author on reasonable request, subject to its permission. With the data, `python run_all.py` reproduces every result (≈ 4 min) with Python 3 and the pinned packages in `requirements.txt`.

# References

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall.

Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and Practice* (3rd ed.). OTexts.

Künsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. *Annals of Statistics*, 17, 1217-1241.

Lawrence, C. (2007). The husbandry of zebrafish (*Danio rerio*): a review. *Aquaculture*, 269, 1-20.

Pedregosa, F., et al. (2011). Scikit-learn: machine learning in Python. *Journal of Machine Learning Research*, 12, 2825-2830.

Schwarz, G. (1978). Estimating the dimension of a model. *Annals of Statistics*, 6, 461-464.

Stumm, W., & Morgan, J. J. (1996). *Aquatic Chemistry* (3rd ed.). Wiley.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society B*, 58, 267-288.

Timmons, M. B., & Ebeling, J. M. (2010). *Recirculating Aquaculture* (2nd ed.). Cayuga Aqua Ventures.

Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261-272.

Zhang, N. R., & Siegmund, D. O. (2007). A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data. *Biometrics*, 63, 22-32.
