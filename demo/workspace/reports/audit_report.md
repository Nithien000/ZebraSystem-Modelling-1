> **SYNTHETIC DEMO DATA - not real facility data.** Results below describe invented data.

# Step 1 - Data audit report

Water-quality export: 316 rows; 0 were deleted in the app (archived) and removed; 0 readings from System 4/5 removed; 0 Aquaneering readings before 8 Jul 2026 removed.
Duplicated days (same system, same date) after that: 0 pairs; the later entry of each pair was kept.

Empty DataFrame
Columns: [system, date, time, temperature, conductivity, ph]
Index: []

## Coverage per system

        system      first       last  calendar_days  days_with_reading  coverage_%  missing_days  longest_gap_days  weekend_readings
   Main System 2026-04-23 2026-11-10            202                128        63.4            74                 6                33
SA Aquaneering 2026-07-08 2026-11-08            124                 68        54.8            56                 5                20
 SA Tecniplast 2026-04-26 2026-11-10            199                120        60.3            79                 5                33

## Reading time of day (the protocol says 09:00)

                count    min    25%    50%    75%    max
system                                                  
Main System     128.0  10.00  10.92  11.63  12.33  12.97
SA Aquaneering   68.0  10.07  10.62  11.43  12.12  12.92
SA Tecniplast   120.0  10.02  10.78  11.33  12.05  12.98

Share of readings between 08:30 and 09:30: 0.0%; after 12:00: 32.6%

## Value ranges and units

               temperature              conductivity                 ph            
                       min median   max          min  median   max  min median  max
system                                                                             
Main System           26.5   28.1  29.3          473   647.0  1058  7.4    7.7  8.0
SA Aquaneering        26.7   28.3  29.3          940  1019.5  1287  7.4    7.7  8.0
SA Tecniplast         26.7   28.2  29.4          530   711.0  1080  7.2    7.6  8.0

Units inferred from the ranges: temperature in °C (26-30), conductivity in µS/cm (560-1570, typical for zebrafish systems), pH in pH units (6.7-9.0). No unit column exists in the export, so this is an inference.

## Flagged readings (possible outliers or real jumps)

       system       date    parameter  value
  Main System 2026-08-20 conductivity 1058.0
SA Tecniplast 2026-04-30  temperature   28.6
SA Tecniplast 2026-10-10 conductivity 1080.0

## Largest day-to-day conductivity changes (between consecutive readings)

        system       date  change_uS  days_since_prev
   Main System 2026-08-22     -492.0                2
   Main System 2026-08-20      487.0                3
   Main System 2026-09-20      260.0                1
   Main System 2026-08-03     -145.0                2
   Main System 2026-06-15      138.0                7
 SA Tecniplast 2026-10-10      446.0                1
 SA Tecniplast 2026-10-11     -391.0                1
 SA Tecniplast 2026-08-21      228.0                2
 SA Tecniplast 2026-10-05     -168.0                2
 SA Tecniplast 2026-06-01       96.0                1
SA Aquaneering 2026-09-01      290.0                3
SA Aquaneering 2026-08-12       74.0                2
SA Aquaneering 2026-08-20      -70.0                1
SA Aquaneering 2026-08-10      -66.0                1
SA Aquaneering 2026-11-08      -62.0                4

## Tank records

Tank export: 28 tanks (0 archived). Fish count per day was rebuilt from each tank's setup date, mortality log and archive date.
Tanks never mapped to one of the three systems (rack name):
Series([], )

Fish per system, monthly mean of the daily total:
system      Main System  SA Tecniplast  SA Aquaneering
date                                                  
2026-04-01          NaN          120.0             NaN
2026-05-01       1280.0          120.0             NaN
2026-06-01       1280.0          120.0             NaN
2026-07-01       1280.0          120.0           800.0
2026-08-01       1280.0          120.0           800.0
2026-09-01       1280.0          239.0           800.0
2026-10-01       1280.0          330.0           800.0
2026-11-01       1280.0          330.0           800.0

Tanks per system, monthly mean:
system      Main System  SA Tecniplast  SA Aquaneering
date                                                  
2026-04-01          NaN            4.0             NaN
2026-05-01          8.0            4.0             NaN
2026-06-01          8.0            4.0             NaN
2026-07-01          8.0            4.0            10.0
2026-08-01          8.0            4.0            10.0
2026-09-01          8.0            7.0            10.0
2026-10-01          8.0           10.0            10.0
2026-11-01          8.0           10.0            10.0

Cleaning events per month (tank cleanings logged in the app):
system   Main System  SA Tecniplast  SA Aquaneering
date                                               
2026-05          8.0            4.0             NaN
2026-06          4.0            4.0             NaN
2026-07          6.0            3.0            10.0
2026-08          8.0            4.0             5.0
2026-09          8.0            3.0             6.0
2026-10          6.0            9.0            10.0
2026-11          1.0            3.0             3.0

Feeding log: present in the app but judged unreliable by the facility manager; NOT used.