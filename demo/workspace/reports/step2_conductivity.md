> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 2 - Conductivity mass-balance model

Fitted on ALL readings of each system, to show the model structure and the inferred
events. The calibration (Apr-Aug) / validation (Sep-Oct) split and the parameter
uncertainty come in Step 4.

## Main System  (volume 525 L)

- Readings used: 127 (2026-04-23 to 2026-11-10)
- Dosing controller ON (inferred from flat stretches): 2026-04-26 to 2026-05-04, 2026-07-02 to 2026-07-15
- Effective renewal rate r = 0.129 per day (facility assumption 0.10); a jump loses half its size in 5.4 days
- Forcing r = 0.10 instead: BIC +7.2 (worse)
- Equilibrium without intervention C* = 638 µS/cm
- Controller set-point C_set = 650 µS/cm (strength k = 0.5 per day)
- Temperature effect on evaporation: BIC +4.8 -> not supported by the data, left out
- Inferred jumps: 3 (2 up = salt additions, 1 down = dilutions)
- In-sample fit: RMSE 24.8 µS/cm, MAE 19.8, R² 0.759
- Check: exact daily solution vs solve_ivp differ by at most 1.0e-07 µS/cm

## SA Tecniplast  (volume 175 L)

- Readings used: 119 (2026-04-26 to 2026-11-10)
- Dosing controller ON (inferred from flat stretches): 2026-06-24 to 2026-07-03, 2026-09-21 to 2026-09-28, 2026-10-17 to 2026-10-29
- Effective renewal rate r = 0.209 per day (facility assumption 0.10); a jump loses half its size in 3.3 days
- Forcing r = 0.10 instead: BIC +31.1 (clearly worse)
- Equilibrium without intervention C* = 709 µS/cm
- Controller set-point C_set = 707 µS/cm (strength k = 0.5 per day)
- Temperature effect on evaporation: BIC +3.7 -> not supported by the data, left out
- Inferred jumps: 2 (1 up = salt additions, 1 down = dilutions)
- In-sample fit: RMSE 22.9 µS/cm, MAE 17.9, R² 0.742
- Check: exact daily solution vs solve_ivp differ by at most 1.1e-07 µS/cm

## SA Aquaneering  (volume 150 L)

- Readings used: 68 (2026-07-08 to 2026-11-08)
- Dosing controller ON (inferred from flat stretches): 2026-07-21 to 2026-08-05, 2026-10-04 to 2026-10-10, 2026-10-27 to 2026-11-01
- Effective renewal rate r = 0.044 per day (facility assumption 0.10); a jump loses half its size in 15.9 days
- Forcing r = 0.10 instead: BIC +26.8 (clearly worse)
- Equilibrium without intervention C* = 986 µS/cm
- Controller set-point C_set = 1026 µS/cm (strength k = 0.5 per day)
- Temperature effect on evaporation: BIC +3.1 -> not supported by the data, left out
- Inferred jumps: 1 (1 up = salt additions, 0 down = dilutions)
- In-sample fit: RMSE 23.9 µS/cm, MAE 18.3, R² 0.913
- Check: exact daily solution vs solve_ivp differ by at most 1.5e-07 µS/cm

## Inferred jumps (NOT observed: estimated from the data)

        system       date  jump_uS_cm                               type  approx_salt_g
   Main System 2026-06-11       207.5           salt addition (inferred)           54.5
   Main System 2026-08-02      -166.0 dilution / water change (inferred)            NaN
   Main System 2026-09-20       238.3           salt addition (inferred)           62.6
 SA Tecniplast 2026-08-20       296.5           salt addition (inferred)           25.9
 SA Tecniplast 2026-10-04      -213.9 dilution / water change (inferred)            NaN
SA Aquaneering 2026-08-30       319.7           salt addition (inferred)           24.0

approx_salt_g: rough conversion 0.5 mg/L per µS/cm (NaCl-like salt) x system volume. An order of magnitude, not a measurement.