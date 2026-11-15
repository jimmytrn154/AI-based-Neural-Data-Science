# P2.1 Handoff: External EEG Microstate Template Matching
Date: 2026-11-15 | Author: Minseo's Colab analysis bundle

## 1. Scope
Primary question: Does our group's K=5 microstate E resemble the scalp topographies of E maps in published template sets?
This is a TOPOGRAPHIC comparison only. It does NOT establish functional equivalence with literature E, cognitive conflict processing, or significance of a hypothesis test.

## 2. What is included
- code/Microstate_E_external_validation_COLAB.ipynb — the executable notebook that generated the external template comparison
- code/colab_pipeline.py — equivalent reproducible Python code
- code/test_parsers.py and code/test_e2e_synthetic.py — offline synthetic tests; supporting fixtures in code/tests/fixtures/ (run the parser test first)
- results/EXTERNAL_MATCHING_ALL_STUDIES.csv — complete 12-reference-set summary, not restricted to positive matches
- results/EXTERNAL_E_SHARED_REST_SIMON_COMPARISON.csv — full vs rest-only vs Simon-only E topographic matches
- results/*_absr.csv — full 5 × K polarity-invariant spatial correlation matrices for each reference set
- results/*_shared_variance_pct.csv — squared spatial correlations for each reference set
- results/*_standard_spline_DIAGNOSTIC.csv — approximate standard 1005 interpolation ONLY, not native-coordinate alignment
- figures/*_heatmap.png — corresponding comparison heatmaps
- original_source_kit/README_EN.md — provenance and implementation notes from the original Colab code kit

Source datasets were downloaded by the notebook from the public MS-Template-Explorer repository; team template data were loaded from the team's public GitHub repository. The original external `.set` files and local copy of `.fif` are omitted from this handoff to keep repository size small. Re-download using the notebook for full reproduction. The separate full Colab output ZIP has copies of the sources.

## 3. Primary method and channel matching
1. Read numeric scalp maps and their channel labels from the source `.set` templates and the team's `.fif` templates.
2. Retain the INTERSECTION of channel *names* for each reference set (19–62 channels across the 12 reference sets). Align vectors in the same channel-name order.
3. For each of our 5 maps against each reference class, compute |Pearson spatial r|. The absolute value disregards global map polarity. Shared variance (%) is 100 * r².
4. Report the closest published class and the published E correlation separately; do not equate them when the closest class is not E.

PRIMARY RESULTS DO NOT USE INTERPOLATION. The `standard_1005_approx_diagnostic` outputs are sensitivity checks using nominal standard montage coordinates; they are not validated source-native electrode coordinate interpolation and must not be represented as such (`native_coordinate_matching_completed=False`).

NO EXPLICIT RE-REFERENCING OPERATION was applied to the map arrays as a separate processing step. Pearson correlation mean-centers the matched channel vectors, so adding/subtracting a common constant reference potential across all matching channels does not affect r. This does NOT solve montage, electrode-location, interpolation, or upstream reference-method differences.

The number of matching electrodes differs across studies; notably the 19-channel and 24-channel comparisons are less directly comparable to dense-montage matches.

## 4. Main findings (our shared K5 map E; |r|)
Zanesco 2020 (K5): best E = 0.8683075, 59 common electrodes; runner-up D = 0.4310257.
Zanesco 2021 meditation (K6, source file `Zanesco 2021a.set`): best E = 0.8432092, 62 channels; runner-up F = 0.7747844.
Hwang 2025 K4 (labels A/B/C/E): best E = 0.5129551, 55 channels.
Custo 2017 K7: best D = 0.7197695, while published E = 0.0798405, 59 channels.
Qin 2022 depression control K5: best C = 0.6422005, published E = 0.3403214, 57 channels.
Qin 2022 depression patient K5: best C = 0.4955447, published E = 0.0988297, 57 channels.
Six additional K4 references lack a published E label and match mostly to C; see complete summary CSV.

CRITICAL literature identification: `Zanesco 2021a.set` is from a meditation-related paper, NOT Zanesco/Denkova/Jha's distinct sustained-attention SART paper. No spatial match to the SART paper's original E map has yet been established.

## 5. Follow-up work required by current team P2.1 plan
- Control matrix for our A–D plus all reference classes (full 5×K matrices are available in outputs).
- Compare to nonmatching-reference baseline correlations / empirically calibrate similarity thresholds.
- Validate native electrode coordinates or jointly defined re-referencing/montage and interpolation strategy; explicitly rerun primary results if necessary.
- Add and evaluate relevant Koenig meta-microstate normative atlas if appropriate.
- Where subject-level maps are available, test stability / bootstrap and consider independent held-out sensitivity checks.

Do not represent all 12 dataset entries as 12 independent cohorts: Qin has two subsets (patients and controls), and source-file labels represent the published template sets, not necessarily independent functional replications.
