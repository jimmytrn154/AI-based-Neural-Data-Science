# P2.1: state E vs. published normative templates

Owner: **Minseo** (analysis, Nov 10–15). Repo integration and verification: Chuong, Nov 15.

Does our K = 5 state E look like the "E" maps in published microstate template sets? This is a **topographic** comparison only: it says nothing about shared generators or function.

## Contents

| Path | What it is |
|---|---|
| `README_EN.md` | Minseo's handoff note (method, findings, caveats) |
| `original_source_kit/README_EN.md` | Notes from the original Colab session |
| `code/Microstate_E_external_validation_COLAB.ipynb` | Minseo's Colab notebook (downloads all inputs itself) |
| `code/colab_pipeline.py` | The same analysis as a Python module; `test_parsers.py`, `test_e2e_synthetic.py` and `tests/fixtures/` are its offline tests (the fixtures hold *our* maps rounded to 1e-4, not reference maps) |
| `results/`, `figures/` | Minseo's outputs: per-dataset \|r\| matrices (`*_absr.csv`), shared variance, interpolation diagnostics, heatmaps, summaries |
| **`p2_1_normative_match.py`** | **Repo version** (Chuong): same method, reads our templates through `src/` and `config.yaml`, adds the checks below; writes `results/features/p2_1_normative_*.csv` and `results/figures/minseo_p2_1_normative_match.png` |

Run from the repo root (the `.set` files are not in the repo; `--download` fetches them):

```bash
python P2_1/p2_1_normative_match.py --refs <folder for the .set files> --download
```

Needs `mat73` (in `requirements.txt`) for the MATLAB v7.3 file (Custo 2017).

## The 12 reference template sets

All from the **MS-Template-Explorer** database (Thomas Koenig, <https://github.com/ThomasKoenigBern/MS-Template-Explorer>, folder `Template Maps/`; repository licence GPL-3.0). Each `.set` file carries its own citation (`EEG.msinfo.Citation`), which is the source of the table below.

**Not committed:** the `.set` files are third-party data, and they also contain the authors' contact e-mails. Download them from the repository above (or run the script with `--download`).

| File | Study | Journal, year, DOI | N | Eyes | Band-pass | Classes |
|---|---|---|---|---|---|---|
| `Custo 2017.set` | Custo et al., Electroencephalographic resting-state networks: source localization of microstates | Brain Connectivity 2017, 10.1089/brain.2016.0476 | 164 | closed | 1–40 Hz | A–G (7) |
| `Zanesco 2020.set` | Zanesco et al., Within and between-person correlates of the temporal dynamics of resting EEG microstates | NeuroImage 2020, 10.1016/j.neuroimage.2020.116631 | 187 | both | 1–45 Hz | A–E (5) |
| `Zanesco 2021a.set` | Zanesco et al., Meditation training modulates brain electric microstates and felt states of awareness | Human Brain Mapping 2021, 10.1002/hbm.25430 | 60 | closed | 40 Hz | A–F (6) |
| `Hwang 2025.set` | Hwang et al., EEG microstate and functional connectivity analyses for differentiating suicide attempt from suicidal ideation in MDD | J Psychiatry Neurosci 2025, 10.1139/jpn-25-0078 | 130 | open | 1–40 Hz | A, B, C, E (4) |
| `Qin2022 Depression_Controls.set` | Qin et al., EEG microstate temporal dynamics predict depressive symptoms in college students (controls) | Brain Topography 2022, 10.1007/s10548-022-00905-0 | 68 | closed | 0.5–45 Hz | A–E (5) |
| `Qin2022 Depression_Patients.set` | Same study (patients) | same | 68 | closed | 0.5–45 Hz | A–E (5) |
| `Zelenina_2022 Oxytocin.set` | Zelenina et al., Temporal dynamics of intranasal oxytocin in human brain electrophysiology | Cerebral Cortex 2022, 10.1093/cercor/bhab404 | 19 | both | 2–20 Hz | A–D (4) |
| `Giordano 2018.set` | Giordano et al., Neurophysiological correlates of avolition-apathy in schizophrenia: a resting-EEG microstates study | NeuroImage: Clinical 2018, 10.1016/j.nicl.2018.08.031 | 206 | closed | 2–20 Hz | A–D (4) |
| `Lassi2023.set` | Lassi et al., Degradation of EEG microstates patterns in subjective cognitive decline and mild cognitive impairment | NeuroImage: Clinical 2023, 10.1016/j.nicl.2023.103407 | 103 | closed | high-pass 1 Hz | A–D (4) |
| `Milz_2016_Cognitive Style.set` | Milz et al., The functional significance of EEG microstates: associations with modalities of thinking | NeuroImage 2016, 10.1016/j.neuroimage.2015.08.023 | 61 | closed | 2–20 Hz | A–D (4) |
| `Schiller2021 Alcohol Intoxication.set` | Schiller et al., Acute alcohol intoxication modulates the temporal dynamics of resting EEG networks | Addiction Biology 2021, 10.1111/adb.13034 | 24 | closed | 2–20 Hz | A–D (4) |
| `Toplutas 2023 DOC_HC.set` | Toplutaş et al., EEG microstate analysis in patients with disorders of consciousness and its clinical significance | Brain Topography 2023, 10.1007/s10548-023-00939-y | 28 | closed | 2–20 Hz | A–D (4) |

The 12 sets come from 11 studies (Qin 2022 contributes two). `Zanesco 2021a` is the meditation study, not Zanesco and colleagues' sustained-attention paper.

## Method (verified Nov 15)

- **Correlation:** centred Pearson, absolute value (polarity ignored). Each map is mean-centred over the channels used, so this is not an uncentred correlation.
- **Electrodes:** matched by channel **name** (case-insensitive); only the channels both montages share are used (at least 15).
- **Re-referencing:** no separate step, but none is needed: centring both maps over the common channels is exactly re-referencing them to the average of that common subset.
- **Interpolation:** published maps splined onto our positions with nominal `standard_1005` coordinates, not each study's native coordinates. **Diagnostic only**; skipped for sets with fewer than 20 usable positions (Toplutas).

**Reproduction:** from the original `.set` files, the repo script reproduces every one of Minseo's 12 primary matrices and all interpolation diagnostics (largest difference 5e-7), including Zanesco 2020 \|r\| = 0.868 and Zanesco 2021a \|r\| = 0.843.

## Results (`results/features/p2_1_normative_*.csv`)

**Common channels with our 64, and the ones missing:**

| Set | Common | Missing from our 64 | Scalp regions missing |
|---|---|---|---|
| Zanesco 2021a | 62 | I1, I2 | inion row |
| Zelenina 2022 | 60 | Iz, I1, I2, TP9 | inion row, left low-temporal |
| Custo 2017 | 59 | Iz, I1, I2, FT10, AFz | inion row, right low-temporal, frontopolar midline |
| Zanesco 2020 | 59 | I1, I2, TP9, TP10, FT10 | inion row, low-temporal |
| Milz 2016 | 58 | Iz, I1, I2, TP9, TP10, FT10 | inion row, low-temporal |
| Qin 2022 (both) | 57 | Iz, I1, I2, TP9, TP10, FT10, AFz | inion row, low-temporal, frontopolar midline |
| Hwang 2025 | 55 | Iz, I1, I2, TP9, TP10, FT10, AF7, AFz, AF8 | inion row, low-temporal, anterior-frontal row |
| Lassi 2023 | 54 | Iz, I1, I2, TP9, TP10, FT10, T7, T8, P7, P8 | inion row, temporal, lateral parietal |
| Schiller 2021 | 51 | F3, F4, F7, F8, C1, C2, P5, P6, PO7, PO8, POz, I1, I2 | lateral frontal, parieto-occipital, inion |
| Giordano 2018 | 24 | 40 channels | most of the 10-10 positions (central, frontal, temporal, parietal) |
| Toplutas 2023 | 19 | 45 channels | most of the 10-10 positions (only a 10-20-like subset left) |

(Full channel lists and region counts per set: `p2_1_normative_E_summary.csv`.)

**Our E against each set** (\|r\|; full 5 × K matrices in `p2_1_normative_matrix.csv`):

| Set | Best match | \|r\| | Second | \|r\| | Gap | Their "E"? | Interpolation diagnostic (best) |
|---|---|---|---|---|---|---|---|
| Zanesco 2020 | E | 0.868 | D | 0.431 | 0.437 | yes | E 0.867 |
| Zanesco 2021a | E | 0.843 | F | 0.775 | 0.068 | yes | E 0.845 |
| Custo 2017 | D | 0.720 | G | 0.670 | 0.050 | no (their E: 0.080) | D 0.720 |
| Zelenina 2022 | C | 0.681 | A | 0.293 | 0.389 | no E class | C 0.649 |
| Qin 2022 controls | C | 0.642 | E | 0.340 | 0.302 | no (their E: 0.340) | C 0.627 |
| Giordano 2018 | C | 0.601 | B | 0.231 | 0.370 | no E class | C 0.576 |
| Schiller 2021 | C | 0.537 | D | 0.364 | 0.174 | no E class | C 0.477 |
| Toplutas 2023 | C | 0.522 | B | 0.148 | 0.374 | no E class | — |
| Hwang 2025 | E | 0.513 | B | 0.275 | 0.238 | yes | E 0.469 |
| Milz 2016 | C | 0.499 | A | 0.191 | 0.308 | no E class | C 0.445 |
| Qin 2022 patients | C | 0.496 | A | 0.221 | 0.275 | no (their E: 0.099) | C 0.454 |
| Lassi 2023 | C | 0.451 | B | 0.228 | 0.222 | no E class | C 0.304 |

- **Count:** our E's best match is the set's E-type map in **3 of 12** sets (Zanesco 2020, Zanesco 2021a, Hwang 2025). Six sets have an E class at all, so 3 of 6. In the six 4-class sets without E, it matches their C.
- **Only the two Zanesco sets reach \|r\| ≥ 0.8**; in Zanesco 2021a the runner-up (F, 0.775) is close.
- **Rest-only and Simon-only E** (`p2_1_normative_rest_simon.csv`) give the same best matches; Zanesco 2020 E 0.815 (Rest-only) and 0.905 (Simon-only), Zanesco 2021a E 0.834 and 0.892.
- **Controls** (`p2_1_normative_controls_AD.csv`): our A and B best-match the same letter in 12 of 12 sets, C in 11, D in 10, at \|r\| 0.83–0.99 (median about 0.95).
- **Baseline** (`p2_1_normative_baseline.csv`): \|r\| between **different** classes within the same published set has a median of 0.42–0.77 and a maximum of 0.68–0.88 (our own C–D: 0.80). An \|r\| around 0.85 is therefore only just above what two distinct classes can reach.

Interpretation is deliberately left out here (see `PROGRESS.md`, P2.1 and P2.5).

## Remaining gaps

- **Koenig et al. 2024 meta-maps** (MetaMaps2023, *Brain Topography*, 10.1007/s10548-023-00993-6) are **not** among the 12; the consolidated atlas has to be obtained from the official MICROSTATELAB/MS-Template-Explorer tools.
- **Custo 2017 (A–G)** is among the 12.
- **Native electrode coordinates** were not verified for any set; a check against the official MATLAB MS-Template-Explorer implementation is still open.
- **Low-density sets:** Giordano (24 channels) and Toplutas (19) share few channels with our montage and are less comparable.
- **Conditions differ:** almost all references are eyes-closed rest; ours are eyes-open rest plus a task; band-passes range from 0.5–45 Hz to 2–20 Hz.
- **Uncertainty:** no confidence interval for \|r\| yet (a bootstrap over our subjects' individual maps would need those maps saved).
- **Labels are study-specific:** a matching letter is not a matching topography (Custo's E is a different map from Zanesco's E).
