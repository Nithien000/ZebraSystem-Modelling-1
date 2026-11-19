> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 4 - Calibration (Apr-Aug), validation (Sep-Oct) and parameter uncertainty

Calibration: readings up to 2026-08-31. Validation: 2026-09-01 to 2026-10-31. Bootstrap: 300 moving-block resamples (block = 5 readings); intervals are 95 % percentile intervals, conditional on the inferred jump dates.

## conductivity

### Main System

- Calibration: 81 readings (2026-04-23 to 2026-08-30); in-sample RMSE 24 µS/cm, R² 0.71; 2 inferred jumps
- Parameters (point estimate [95 % CI]):
    - r: 0.1291 [0.0981, 0.1779]
    - half_life_days: 5.369 [3.896, 7.066]
    - Cstar: 635.1 [628.4, 642.5]
    - Cset: 648.9 [639.3, 664.2]
- Validation, next-reading forecast (38 forecasts): model RMSE 56.6, MAE 38, R² 0.07 | persistence RMSE 62.9, MAE 44.6, R² -0.14 -> model better
- Stability: r refitted on Sep-Oct only = 0.164 per day; calibration CI [0.0981, 0.178] -> inside (stable)

### SA Tecniplast

- Calibration: 72 readings (2026-04-26 to 2026-08-27); in-sample RMSE 24.3 µS/cm, R² 0.71; 1 inferred jumps
- Parameters (point estimate [95 % CI]):
    - r: 0.1853 [0.1245, 0.2543]
    - half_life_days: 3.741 [2.726, 5.567]
    - Cstar: 705.5 [699.4, 711.5]
    - Cset: 701.3 [679.2, 724.4]
- Validation, next-reading forecast (41 forecasts): model RMSE 35.1, MAE 21, R² 0.40 | persistence RMSE 40.2, MAE 26.4, R² 0.21 -> model better
- Stability: r refitted on Sep-Oct only = 0.266 per day; calibration CI [0.124, 0.254] -> outside (parameter changed or poorly determined)

### SA Aquaneering

- Calibration: 30 readings (2026-07-08 to 2026-08-29); in-sample RMSE 20.5 µS/cm, R² 0.07; 0 inferred jumps
- Parameters (point estimate [95 % CI]):
    - r: 1 [0.005, 1]
    - half_life_days: 0.6931 [0.6931, 138.6]
    - Cstar: 998.8 [830.3, 1001]
    - Cset: 1043 [1006, 1068]
- Validation, next-reading forecast (33 forecasts): model RMSE 108, MAE 82.3, R² -0.58 | persistence RMSE 58.6, MAE 33, R² 0.53 -> persistence better
- Stability: r refitted on Sep-Oct only = 0.0387 per day; calibration CI [0.005, 1] -> inside (stable)

## pH

### Main System

- Calibration: 63 readings (2026-05-22 to 2026-08-30); in-sample RMSE 0.112 pH, R² 0.12; 0 inferred jumps; terms: relaxation only
- Parameters (point estimate [95 % CI]):
    - r: 0.01 [0.005, 0.2013]
    - half_life_days: 69.31 [3.455, 138.6]
    - pH_eq_or_in: 7.837 [7.69, 8.114]
- Validation, next-reading forecast (38 forecasts): model RMSE 0.128, MAE 0.094, R² 0.09 | persistence RMSE 0.13, MAE 0.0947, R² 0.07 -> about equal
- Stability: r refitted on Sep-Oct only = 0.204 per day; calibration CI [0.005, 0.201] -> outside (parameter changed or poorly determined)

### SA Tecniplast

- Calibration: 72 readings (2026-04-26 to 2026-08-27); in-sample RMSE 0.144 pH, R² 0.02; 0 inferred jumps; terms: relaxation only
- Parameters (point estimate [95 % CI]):
    - r: 0.5298 [0.005, 1]
    - half_life_days: 1.308 [0.6931, 138.6]
    - pH_eq_or_in: 7.59 [7.397, 7.704]
- Validation, next-reading forecast (42 forecasts): model RMSE 0.142, MAE 0.112, R² -0.61 | persistence RMSE 0.165, MAE 0.138, R² -1.18 -> model better
- Stability: r refitted on Sep-Oct only = 1 per day; calibration CI [0.005, 1] -> inside (stable)

### SA Aquaneering

- Calibration: 30 readings (2026-07-08 to 2026-08-29); in-sample RMSE 0.103 pH, R² 0.35; 1 inferred jumps; terms: relaxation only
- Parameters (point estimate [95 % CI]):
    - r: 0.01 [0.005, 0.1178]
    - half_life_days: 69.31 [5.903, 138.6]
    - pH_eq_or_in: 7.168 [6.222, 7.663]
- Validation, next-reading forecast (33 forecasts): model RMSE 0.146, MAE 0.117, R² -0.84 | persistence RMSE 0.149, MAE 0.118, R² -0.91 -> about equal
- Stability: r refitted on Sep-Oct only = 0.01 per day; calibration CI [0.005, 0.118] -> inside (stable)
