# Microstate E External Topographic Matching — Original Analysis Notes and Reproduction Kit

Original notes dated 2026-11-10 (Asia/Seoul). This document is an English translation of the original source-kit README; reported execution status refers to that original session. For the latest 12-reference handoff and interpretation, see `../README_EN.md`.

## Before you start
**No additional files are required for the initial matching workflow.** The team's original `results/templates/templates_K5.fif` on GitHub includes the names and 3-D positions of the 64 EEG channels. The Colab notebook downloads this file directly.

**Reproduce in Colab:** upload `Microstate_E_external_validation_COLAB.ipynb` to Google Colab, select **Runtime > Run all**, and retrieve the output files. The user must authenticate and run Colab in their own browser; the original analysis session could not operate those browser controls on the user's behalf.

## Work verified in the original session

| Task | Status | Result |
|---|---|---|
| Verify E vs. A–D correlations in an earlier user ZIP | Calculated | E–C `|r|=0.336651` |
| Reproduce K4 ↔ K5 template comparison | Calculated | K4 Central–K5 E `|r|=0.985921` |
| Reproduce Rest/Simon template comparison | Calculated | Rest-only E–Simon-only E `|r|=0.948385` |
| Read the actual K5 `.fif` from the team GitHub repository | Executed | Extracted all 64 double-precision E-map values and 64 three-dimensional channel positions |
| Verify exact K5 E against rounded values in the earlier ZIP | Calculated | Maximum difference per channel: 4.704e-05 |
| Compare our exact E with the preserved published Hwang 2025 E map | Calculated | 55 shared channels; signed r ≈ −0.513004; `|r|=0.513004`; shared variance 26.32% |
| Compare our exact E with the preserved published Qin 2022 control C map | Calculated | 57 shared channels; `|r|=0.642195`; shared variance 41.24% |
| Approximate spline sensitivity check with MNE `standard_1005` | Calculated | Hwang E `|r|=0.467162`; Qin C `|r|=0.623439`; **not** native electrode-coordinate matching |
| Implement FIF/EEGLAB `.set` parsers and test synthetic inputs | Code tested | Verified handling of K5 labels, channels, values, and published Hwang A/B/C/E labels |
| Implement Colab workflow to download templates, compute full matrices, compare labels, produce plots/CSVs | Code written; syntax and synthetic tests run | **The original session had not yet run the full internet-connected Colab workflow or parsed the Custo source.** The subsequently supplied `results/` directory contains later execution outputs; consult those separately. |

## Analyses needing additional data

- **Participant bootstrap:** individual K5 A–E maps were not present under `results/templates/` on the team GitHub repository. Request `subject_maps.npz` with `maps` of shape `(n_subjects, 5, 64)` and optional `groups` of shape `(n_subjects,)`. The code uses Hungarian matching to align each participant's maps. A bootstrap of stored maps is approximate, not a bootstrap that refits the original two-level modified k-means from raw EEG.
- **Participant-level Rest/Simon TANOVA:** requires per-participant, per-condition maps. Proposed input `subject_maps_rest_simon.npz`: `rest_maps` and `simon_maps`, each `(n_subjects, 5, 64)` in identical subject order, plus `subject_ids` and `groups`. A group-average published reference map alone does not support a two-sample cross-cohort equivalence test; individual data from the reference cohort would be needed.
- **External-template back-fitting to EEG and held-out GEV:** requires preprocessed raw EEG or participant-level GFP-peak channel vectors and preprocessing/channel/reference metadata. Public team feature CSVs contain summaries, not the raw input needed to refit or back-fit these maps.
- **Rigorous native-coordinate cross-study matching:** requires verifying actual electrode coordinates, units, coordinate frames, and spherical interpolation conventions in each published `.set` file. The Colab script reports approximate sensitivity results using nominal `standard_1005` coordinates, and explicitly does not claim complete native-coordinate alignment. The results should be checked against the official MATLAB MS-Template-Explorer implementation.

## Source-kit file inventory (some files are not part of this smaller handoff)

- `Microstate_E_external_validation_COLAB.ipynb`: self-contained notebook that downloads FIF templates from the team repository and multiple `.set` files from the public reference repository.
- `colab_pipeline.py`: reusable Python implementation of the notebook's analysis functions.
- `verify_offline.py`: independent numerical verification for two preserved external reference vectors from the original session (available in the original reproduction kit, not in this reduced handoff).
- `results/internal_K5_AtoE_absr.csv`, `results/K4_vs_K5_absr.csv`: original internal-validation outputs (in the original reproduction kit).
- `results/external_map_checks_exact_our_E.csv`: original partial checks against Hwang and Qin reference maps.
- `inputs/k5_chan_coords_from_original_fif.csv`, `inputs/our_K5_E_exact_from_github_fif.csv`: values extracted from the team's FIF in the original session.
- `figures/Hwang_E_scatter.png`, `figures/Qin_HC_C_scatter.png`: two original comparison figures.
- `test_parsers.py`: tests of the FIF and `.set` parsers using synthetic source files.

## Reproducing in Google Colab

1. Open Colab and upload `Microstate_E_external_validation_COLAB.ipynb`.
2. Select **Runtime > Run all** to install libraries and fetch source data.
3. Inspect `EXTERNAL_MATCHING_ALL_STUDIES.csv` and **also** `REFERENCES_NOT_PARSED.csv`. Never claim a paper was tested when its template could not be parsed.
4. Inspect `EXTERNAL_E_SHARED_REST_SIMON_COMPARISON.csv` to test whether matching patterns recur with the shared, Rest-only, and Simon-only E maps.
5. Inspect `*_shared_variance_pct.csv` and `*_heatmap.png` to compare all of our A–E maps against reference classes A–G.
6. Interpret `native_coordinate_matching_completed=False` correctly: this is **not** confirmation of an official native-coordinate-aligned MS-Template-Explorer comparison.
7. After verifying the results, retrieve `E_Microstate_External_Test.zip`.

## Defensible current conclusion

> Our state E is clearly distinguishable from A–D within the current K=5 solution and reproducible across Rest and Simon. Preliminary cross-study comparisons show that the most similar published map can carry either label C or label E. Whether our E corresponds to the normative E class remains unresolved pending full published-map matching and validated electrode-position alignment.

Matching alphabetic labels ≠ matching scalp topographies ≠ matching neural networks.

## Sources and public reference material

- MICROSTATELAB: Kalburgi et al., *Brain Topography* (2024), DOI [10.1007/s10548-023-01003-5](https://doi.org/10.1007/s10548-023-01003-5), especially the paper's discussion of map sorting, comparison, Fig. 8, and TANOVA.
- Koenig et al., *EEG-Meta-Microstates* (2024), DOI [10.1007/s10548-023-00993-6](https://doi.org/10.1007/s10548-023-00993-6).
- [Team template repository](https://github.com/jimmytrn154/AI-based-Neural-Data-Science/tree/main/results/templates).
- [Published reference templates](https://github.com/ThomasKoenigBern/MS-Template-Explorer/tree/main/Template%20Maps).

**Important note about MetaMaps2023:** The Colab notebook lists public templates from individual published studies, including Custo, not the consolidated MetaMaps2023 atlas itself. If the integrated MetaMaps2023 map vectors and electrodes are needed, they must be obtained from the official program; the original kit's `REQUEST_ADDITIONAL_DATA_KO.md` (not included in this reduced handoff) identified this as an outstanding request.
