> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 5 - Transfer test on SA Aquaneering

(A) full transfer of donor process parameters, (B) donor rate only, (C) Aquaneering's own fit.
System-specific jumps are always re-inferred on Aquaneering.

## conductivity

- (C) Aquaneering own fit: r = 0.0437, C* = 986; RMSE 23.9 µS/cm, BIC 461.1, 1 jumps; forecast RMSE 108 (persistence 58.6)
- (A) from Main System: r = 0.129 (transferred), equilibrium = 638.3 (transferred); RMSE 333 µS/cm, BIC 806.7 (Δ vs C +345.6), 1 jumps; forecast RMSE 122
- (B) from Main System: r = 0.129 (transferred), equilibrium = 1016 (Aquaneering's own); RMSE 24.9 µS/cm, BIC 483.9 (Δ vs C +22.8), 3 jumps; forecast RMSE 57.9
- (A) from SA Tecniplast: r = 0.209 (transferred), equilibrium = 708.7 (transferred); RMSE 344 µS/cm, BIC 798.4 (Δ vs C +337.3), 0 jumps; forecast RMSE 134
- (B) from SA Tecniplast: r = 0.209 (transferred), equilibrium = 1017 (Aquaneering's own); RMSE 26.6 µS/cm, BIC 517.8 (Δ vs C +56.7), 5 jumps; forecast RMSE 60.7

## pH

- (C) Aquaneering own fit: r = 0.149, pH_eq = 7.71; RMSE 0.12 pH, BIC -275.5, 0 jumps; forecast RMSE 0.146 (persistence 0.149)
- (A) from Main System: r = 0.01 (transferred), equilibrium = 7.783 (transferred); RMSE 0.124 pH, BIC -279.6 (Δ vs C -4.0), 0 jumps; forecast RMSE 0.148
- (B) from Main System: r = 0.01 (transferred), equilibrium = 7.701 (Aquaneering's own); RMSE 0.121 pH, BIC -278.7 (Δ vs C -3.2), 0 jumps; forecast RMSE 0.146
- (A) from SA Tecniplast: r = 0.174 (transferred), equilibrium = 7.542 (transferred); RMSE 0.159 pH, BIC -220.5 (Δ vs C +55.0), 2 jumps; forecast RMSE 0.135
- (B) from SA Tecniplast: r = 0.174 (transferred), equilibrium = 7.712 (Aquaneering's own); RMSE 0.12 pH, BIC -279.7 (Δ vs C -4.2), 0 jumps; forecast RMSE 0.116

## Summary table

    variable         donor         variant     r  equilibrium  fit_RMSE      BIC  dBIC_vs_C  jumps  forecast_RMSE  persistence_RMSE
conductivity             -         C refit 0.044      986.046    23.886  461.102      0.000      1        107.867            58.638
conductivity   Main System A full transfer 0.129      638.346   332.878  806.736    345.633      1        122.115            58.638
conductivity   Main System     B rate only 0.129     1015.522    24.949  483.903     22.800      3         57.902            58.638
conductivity SA Tecniplast A full transfer 0.209      708.747   343.530  798.361    337.258      0        133.986            58.638
conductivity SA Tecniplast     B rate only 0.209     1017.499    26.571  517.788     56.686      5         60.653            58.638
          pH             -         C refit 0.149        7.712     0.120 -275.532      0.000      0          0.146             0.149
          pH   Main System A full transfer 0.010        7.783     0.124 -279.556     -4.024      0          0.148             0.149
          pH   Main System     B rate only 0.010        7.701     0.121 -278.723     -3.191      0          0.146             0.149
          pH SA Tecniplast A full transfer 0.174        7.542     0.159 -220.536     54.996      2          0.135             0.149
          pH SA Tecniplast     B rate only 0.174        7.712     0.120 -279.734     -4.201      0          0.116             0.149