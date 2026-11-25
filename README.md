# KEPLER-X v0.5 — Geometry and CVaR benchmark

Exploratory extension of v0.4, after observing its results. Research prototype, not flight navigation.

## Run on Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe geometry_cvar.py
```

To reproduce the bundled 2028 case without network access:

```powershell
.\.venv\Scripts\python.exe geometry_cvar.py --years 2028 --offline
```

The other windows require an Internet connection or cached JPL Horizons files with the exact filenames used by the v0.4 program. Copy the previous `horizons_raw` JSON files into this version's `horizons_raw` folder to reuse them. The bundled 2028 inputs are archived unchanged. Outputs are under `results_v05` (overwritten for each run; preserve copies before reruns).

## Method

Search the same broad launch intervals as v0.4, using 3-day departure and 5-day flight-time spacing, explicitly including 321 days, and refine around the three best coarse candidates under the C3 limit (or lowest-cost points if none are under it). Report the minimum finite delta-v regardless of feasibility and the minimum under the C3 cap. This does not establish global infeasibility.

Compare fixed arrival and arrival flexibility of +/- 7 days relative to the refined nominal arrival, each with its own nominal baseline. Search departures from 28 days before to 14 days after the reference departure, within the launch search interval. Nominal transfer durations and delayed durations must be 140–321 days. Model uses C3 <= 25 km^2/s^2 and total ideal impulsive delta-v <= 6.2 km/s. CVaR90 policies must be feasible for every design delay 0–14 days, with nominal delta-v within 0.1 km/s of the respective baseline. The CVaR is the mean of the worst 10% probability mass, including fractional boundary mass for equally weighted discrete scenarios. Fixed buffers are included when nominally feasible and available in the search domain, even if outside the CVaR penalty cap. Report missing policies rather than silently replacing them.

Test delays 15–30 days are evaluated after policy selection. They are **not independent confirmatory data** because v0.4 results already informed the revised design. Delay distributions and budget are illustrative, not operational probabilities. Geometry outputs give the unsigned angle between planetary position vectors and inclination of the computed prograde transfer plane; no new geometric penalty is yet used in the optimizer. Missing or invalid test outcomes count as infeasible. The underlying Lambert solver is unchanged; an independent solver benchmark remains pending. No multi-body forces, multirevolution solutions, finite burns, real vehicle constraints or navigation validation are included.

## Outputs

- `summary_v05.json`: search diagnostics, source hashes and policies.
- `candidates.csv`: complete nominally feasible design profiles and CVaR metrics.
- `policies.csv`: policy evaluation metrics.
- `geometry.csv`: delay-level C3, delta-v, feasibility and transfer geometry.

## Local test

2028 ran successfully using the supplied Horizons cache: 2,356 coarse evaluations, 201 fine evaluations, 539 nominally feasible candidates with complete design profiles. Fixed CVaR equals the seven-day buffer; flexible CVaR moves arrival one day later and meets 25/31 tested delays vs 24/31 for fixed CVaR. Fourteen-day buffer meets 31/31 but exceeds the illustrative 0.1 km/s nominal penalty cap. Later windows are not yet executed here. Do not interpret this as evidence of general superiority or novelty.

## Repository publication scope

This repository contains the unchanged v0.5 Python source, the archived 2028 Horizons inputs, and the four-window v0.4 and v0.5 summary reports. The complete downloadable publication bundle also contains detailed 2028 benchmark CSVs and earlier v0.4 plots; those supplementary files are not included in this initial GitHub commit. The 2028 detailed outputs can be regenerated offline using the instructions above. The research report remains exploratory and is not peer reviewed. No NASA affiliation or endorsement is claimed. No license has been selected.
