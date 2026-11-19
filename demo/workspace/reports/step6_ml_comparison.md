> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 6 - ML vs mechanistic vs persistence (next-reading forecasts, test Sep-Oct 2026)

RMSE on the test period; in brackets the 95 % bootstrap interval of (method RMSE - persistence RMSE).
Negative = better than persistence. An interval containing 0 = not distinguishable from persistence.

## conductivity

### Main System  (train 78 rows, test 38 rows)

- Persistence: RMSE 62.9
- Ridge (per system): RMSE 58.2  [-12.5, -0.0852] -> better
- Gradient boosting (per system): RMSE 55.4  [-18.5, -1.13] -> better
- Ridge (pooled): RMSE 57.7  [-17.1, +1.07] -> no clear difference
- Gradient boosting (pooled): RMSE 58.5  [-13.4, +1.49] -> no clear difference
- Mechanistic ODE (Step 4): RMSE 56.6  [-15.1, -1.45] -> better

### SA Tecniplast  (train 69 rows, test 41 rows)

- Persistence: RMSE 40.2
- Ridge (per system): RMSE 39.5  [-8.68, +6.21] -> no clear difference
- Gradient boosting (per system): RMSE 43  [-4.33, +10.1] -> no clear difference
- Ridge (pooled): RMSE 40.2  [-4.2, +3.19] -> no clear difference
- Gradient boosting (pooled): RMSE 35.5  [-10.9, +0.227] -> no clear difference
- Mechanistic ODE (Step 4): RMSE 35.1  [-12, -0.987] -> better

### SA Aquaneering  (train 27 rows, test 33 rows)

- Persistence: RMSE 58.6
- Ridge (per system): RMSE 203  [+93.9, +199] -> worse
- Gradient boosting (per system): RMSE 69.9  [-1.18, +31.7] -> no clear difference
- Ridge (pooled): RMSE 58.9  [-3.72, +4.07] -> no clear difference
- Gradient boosting (pooled): RMSE 61.2  [-2.62, +12] -> no clear difference
- Mechanistic ODE (Step 4): RMSE 108  [+26.8, +80.3] -> worse

## pH

### Main System  (train 79 rows, test 38 rows)

- Persistence: RMSE 0.13
- Ridge (per system): RMSE 0.143  [-0.0236, +0.0561] -> no clear difference
- Gradient boosting (per system): RMSE 0.159  [-0.0121, +0.0756] -> no clear difference
- Ridge (pooled): RMSE 0.127  [-0.0334, +0.0306] -> no clear difference
- Gradient boosting (pooled): RMSE 0.14  [-0.0273, +0.05] -> no clear difference
- Mechanistic ODE (Step 4): RMSE 0.128  [-0.00267, -0.000184] -> better

### SA Tecniplast  (train 69 rows, test 42 rows)

- Persistence: RMSE 0.165
- Ridge (per system): RMSE 0.138  [-0.0618, +0.00126] -> no clear difference
- Gradient boosting (per system): RMSE 0.162  [-0.0384, +0.034] -> no clear difference
- Ridge (pooled): RMSE 0.141  [-0.058, +0.0071] -> no clear difference
- Gradient boosting (pooled): RMSE 0.164  [-0.0325, +0.0322] -> no clear difference
- Mechanistic ODE (Step 4): RMSE 0.142  [-0.0452, -0.00422] -> better

### SA Aquaneering  (train 27 rows, test 33 rows)

- Persistence: RMSE 0.149
- Ridge (per system): RMSE 0.126  [-0.0396, -0.00558] -> better
- Gradient boosting (per system): RMSE 0.125  [-0.0563, +0.00851] -> no clear difference
- Ridge (pooled): RMSE 0.111  [-0.0571, -0.0161] -> better
- Gradient boosting (pooled): RMSE 0.117  [-0.0532, -0.0102] -> better
- Mechanistic ODE (Step 4): RMSE 0.146  [-0.00643, +0.000948] -> no clear difference
