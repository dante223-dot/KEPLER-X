# KEPLER-X: Geometric sensitivity and tail-risk analysis of Earth–Mars transfer windows

**Exploratory research report — v0.5, November 2026**

Proposed author: Manuela Maria Romeiro (independent project, Brazil). Authorship and affiliation must be confirmed before release.

## Abstract

This exploratory computational study examines the trade-off between nominal propulsion cost and launch-delay tolerance for Earth–Mars transfers. Geometric heliocentric state vectors from JPL Horizons are used in a prograde zero-revolution Lambert model with ideal impulsive patched-conic departure and capture. Four launch-window search domains were examined, with 10,144 coarse and refined grid evaluations reported. Candidate policies compare nominal cost minimization, fixed launch buffers and conditional value-at-risk (CVaR90) with fixed or flexible arrival dates. Under illustrative constraints, flexible-arrival CVaR increased feasible delay scenarios in the 2028 window from 17/31 (nominal) and 24/31 (seven-day buffer) to 25/31. In 2033 and 2035, feasible scenario counts matched the seven-day buffer. No feasible grid point was identified for the 2031 domain. The results demonstrate trade-offs in a simplified model, not general superiority, operational mission reliability or an established novel optimization method.

## Methods

The model uses Sun-only Lambert transfers, with prograde, zero-revolution arcs and ideal impulsive patched conics. Circular parking orbits of 300 km are assumed. Limits are C3 <= 25 km²/s², total ideal delta-v <= 6.2 km/s and transfer duration 140–321 days. The design delays are equally weighted integer days 0–14; later delays 15–30 are evaluated after policy selection. CVaR90 is the mean of the upper 10% probability mass with fractional weighting at the discrete boundary. CVaR policies must satisfy all design scenarios and limit nominal delta-v penalty to 0.1 km/s relative to their own domain baseline. Flexible arrivals are chosen within ±7 days of the reference date; the selected arrival remains fixed when evaluating delays. The comparison is exploratory because earlier results informed v0.5.

## Reported results

| Window | Policy | Departure | Arrival | Nominal delta-v (km/s) | Design CVaR90 (km/s) | Feasible delays 0–30 |
|---|---|---|---|---:|---:|---:|
| 2028 | fixed_nominal | 2028-11-23 | 2029-09-19 | 5.861254 | 6.047099 | 17/31 |
| 2028 | fixed_buffer_7 | 2028-11-16 | 2029-09-19 | 5.890746 | 5.895783 | 24/31 |
| 2028 | fixed_buffer_14 | 2028-11-09 | 2029-09-19 | 5.971700 | 5.967044 | 31/31 |
| 2028 | fixed_cvar90 | 2028-11-16 | 2029-09-19 | 5.890746 | 5.895783 | 24/31 |
| 2028 | flex_cvar90 | 2028-11-16 | 2029-09-20 | 5.892986 | 5.894238 | 25/31 |
| 2033 | fixed_nominal | 2033-04-18 | 2033-11-04 | 6.051947 | 6.198029 | 14/31 |
| 2033 | fixed_buffer_7 | 2033-04-11 | 2033-11-04 | 6.100519 | 6.097379 | 21/31 |
| 2033 | fixed_cvar90 | 2033-04-11 | 2033-11-04 | 6.100519 | 6.097379 | 21/31 |
| 2033 | flex_cvar90 | 2033-04-11 | 2033-11-03 | 6.094345 | 6.094273 | 21/31 |
| 2035 | fixed_nominal | 2035-06-29 | 2036-01-16 | 5.746832 | 5.859692 | 28/31 |
| 2035 | fixed_buffer_7 | 2035-06-22 | 2036-01-16 | 5.766952 | 5.773310 | 31/31 |
| 2035 | fixed_buffer_14 | 2035-06-15 | 2036-01-16 | 5.834449 | 5.829964 | 31/31 |
| 2035 | fixed_cvar90 | 2035-06-21 | 2036-01-16 | 5.773332 | 5.771579 | 31/31 |
| 2035 | flex_cvar90 | 2035-06-21 | 2036-01-15 | 5.770024 | 5.770332 | 31/31 |

For the 2031 search domain, the minimum scanned delta-v under the C3 cap was 6.217527 km/s, exceeding the illustrative limit by 17.527 m/s. This is not proof of global infeasibility. The 2028 results also show that the nominal transfer plane becomes more inclined across late-delay cases; the reported maximum is 26.744 degrees, versus 7.298 degrees for flexible CVaR. This is a diagnostic association, not a validated causal predictor. No new geometric penalty is used in the optimizer.

## Limitations and interpretation

Fixed CVaR matches a seven-day buffer in 2028 and 2033. A fourteen-day buffer achieves 31/31 in 2028 but lies outside the CVaR nominal penalty cap; this is not an equal-budget comparison. The delay counts combine design and later evaluation scenarios and are not real-world probabilities. The four windows are not independent statistical samples. There is no independent Lambert solver benchmark, multi-body trajectory validation, launch vehicle model, or demonstration of novelty against the literature. The analysis is not suitable for flight navigation. The supplied full-window summary records the user's local execution; only 2028 raw ephemerides are included in the initial GitHub repository; detailed 2028 outputs are available in the separate publication bundle. Later raw response files and detailed outputs are not archived in this bundle; source URLs and hashes are recorded in the summary. Prior-version graphs are labeled v0.4 and should not be presented as v0.5 plots.

## Data and reproducibility

See README.md for execution instructions and the scope of this GitHub publication, and reported_results/summary_v05.json for the four-window run. The complete publication bundle also contains benchmark_2028 and prior-version plots; those supplemental files are not in the initial GitHub commit. The raw JPL input files for 2028 are included. Re-fetching later ephemerides may produce different response metadata or values, so exact byte-level reproducibility of later windows is not asserted. The code is frozen at v0.5; this archive adds reporting material only.

## AI assistance and independence

The code and documentation were developed with substantial generative AI assistance. The human project organizer executed experiments locally and provided outputs and feedback. The report is a draft for human verification and approval; the human author is responsible for any released claims. This project uses publicly available JPL data and is not affiliated with, endorsed by, or validated by NASA/JPL.

## Data source

NASA/JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/
API documentation: https://ssd-api.jpl.nasa.gov/doc/horizons.html

This report is not peer reviewed. A full literature review is outside its current scope.