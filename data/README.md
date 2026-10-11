# Data

This pipeline expects a tab-separated, header-less text file named
`ascii_correct.txt` in the working directory, with exactly nine columns in
this order:

| Column | Field    | Description                                  |
|-------:|----------|-----------------------------------------------|
| 1      | `long`   | Longitude (decimal degrees)                   |
| 2      | `lat`    | Latitude (decimal degrees)                    |
| 3      | `year`   | Origin year                                    |
| 4      | `month`  | Origin month (1-12)                            |
| 5      | `day`    | Origin day of month                            |
| 6      | `mag`    | Magnitude                                      |
| 7      | `sec`    | Origin seconds (may include decimals)          |
| 8      | `hour`   | Origin hour, 0-23 (24-hour clock)              |
| 9      | `minute` | Origin minute, 0-59                            |

## Provenance

The catalog used in the original analysis was obtained from the Iranian
Seismological Center (IRSC). Raw catalog files are **not redistributed in
this repository**; consult IRSC's own data access policy to obtain the
underlying dataset, or substitute your own catalog in the same format.

## Data quality notes

An earlier version of this catalog (used in a prior analysis of the same
case study) had a corrupted `hour` field — every event's timestamp
resolved to a constant hour value, which is not physically plausible for a
12-year catalog and which corrupts the chronological ordering of same-day
events. The pipeline's Step 0 validation (`length(unique(raw$hour)) > 5`)
is a guard against silently re-introducing this issue; if it fails, check
the raw catalog export for a similar timestamp truncation bug (commonly
caused by a spreadsheet cell-format error when exporting `HH:MM:SS` values
as `MM:SS`).
