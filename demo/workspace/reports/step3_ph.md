> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 3 - pH model

Fitted on all readings from the date each system's fish count is reliable. The
calibration/validation split and uncertainty come in Step 4.

## Main System  (volume 525 L, 2.4-2.4 fish per litre)

- Readings used: 109 (2026-05-22 to 2026-11-10)
- Terms kept: renewal only
- BIC with vs without each term (negative = the term helps; kept only if < -2): fish +inf, cleaning +4.5, dosing +7.6
- pH relaxation rate r = 0.010 per day (half-return time 69.3 days). For pH this includes CO2 degassing by aeration, not only water renewal.
- Equilibrium pH (no detectable fish effect) = 7.78
- Other inferred jumps: 0
- In-sample fit: RMSE 0.121, MAE 0.097, R² 0.076
- Check vs solve_ivp: max difference 1.5e-14

## SA Tecniplast  (volume 175 L, 0.7-1.9 fish per litre)

- Readings used: 120 (2026-04-26 to 2026-11-10)
- Terms kept: renewal only
- BIC with vs without each term (negative = the term helps; kept only if < -2): fish -0.2, cleaning +4.4, dosing +6.8
- pH relaxation rate r = 0.174 per day (half-return time 4.0 days). For pH this includes CO2 degassing by aeration, not only water renewal.
- Equilibrium pH (no detectable fish effect) = 7.54
- Other inferred jumps: 1
- In-sample fit: RMSE 0.130, MAE 0.105, R² 0.155
- Check vs solve_ivp: max difference 1.3e-10

## SA Aquaneering  (volume 150 L, 5.3-5.3 fish per litre)

- Readings used: 68 (2026-07-08 to 2026-11-08)
- Terms kept: renewal only
- BIC with vs without each term (negative = the term helps; kept only if < -2): fish +inf, cleaning +0.8, dosing +8.0
- pH relaxation rate r = 0.149 per day (half-return time 4.7 days). For pH this includes CO2 degassing by aeration, not only water renewal.
- Equilibrium pH (no detectable fish effect) = 7.71
- Other inferred jumps: 0
- In-sample fit: RMSE 0.120, MAE 0.094, R² 0.017
- Check vs solve_ivp: max difference 1.1e-10

## Other inferred pH jumps (NOT observed)

       system       date  jump_pH
SA Tecniplast 2026-07-07     0.39