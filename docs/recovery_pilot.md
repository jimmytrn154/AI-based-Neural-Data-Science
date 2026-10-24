# Recovery pilot: go/no-go note (task C4)

Exploratory. Computed Oct 23 by Chuong with `notebooks/chuong_recovery.py`. Per-subject values are in `results/features/recovery.csv` and the figure in `results/figures/chuong_recovery.png`. No group or clinical analyses were run.

## What was measured

140 subjects whose Simon recording passes the recording QC (K = 5 templates).

- **Trials:** responded, non-bad epochs; correct and error trials both need RT 150–1500 ms.
- **Preprocessing:** 2–20 Hz, 250 Hz, then back-fit with the K = 5 templates.
- **Baseline:** each subject's mean state-probability vector (A–E), cue-locked −200 to 0 ms, correct trials.
- **Recovery curve:** distance (Euclidean, between 5-state probability vectors) from that baseline at each time point after the keypress.
- **Recovery measure:** mean distance over +200 to +1000 ms after the keypress.
- **Comparisons, within subject:**
  - incongruent vs. congruent correct trials (needs ≥ 40 per condition; 137 subjects);
  - errors vs. correct trials (needs ≥ 10 errors and at least as many correct as error trials; 111 subjects).
- **Error trials per subject:** median 25 (IQR 15–49, range 1–133).
- **Reliability:** split-half, odd vs. even trial numbers, each half with its own baseline; Spearman–Brown corrected. The config keeps a measure only if reliability ≥ 0.5.

## Result

| Measure (column in `recovery.csv`) | n | Split-half r | Spearman–Brown | Rank-based (Spearman–Brown of ρ) | Decision |
|---|---|---|---|---|---|
| Recovery, all correct trials (`simon_recovery_correct`) | 140 | 0.94 | **0.97** | 0.81 | **Keep** |
| Conflict effect, incongruent − congruent (`simon_recovery_conflict`) | 137 | 0.14 | **0.24** | 0.25 | **Drop** |
| Error effect, error − matched correct (`simon_recovery_error_effect`) | 108 | 0.52 | **0.68** | 0.57 | **Keep (marginal)** |

**Go** for the overall recovery measure and the error effect; **no-go** for the conflict effect.

**Descriptive within-subject differences:**
- incongruent − congruent: median +0.005 (60% of subjects > 0);
- error − matched correct: median +0.020 (67% of subjects > 0).

## Caveats

- **Trial-count bias:** a state-probability vector estimated from fewer trials lies further from the baseline by chance. Errors are therefore compared with random subsets of correct trials of the same size (200 subsets per subject), not with all correct trials. The overall measure does not depend on the number of correct trials (Spearman ρ = 0.03 with trial count).
- **Epoch edge:** the response-locked epochs end at +1000 ms. Back-fitting leaves each epoch's last segment unlabeled, so only 47% of trials are labeled at +900 ms and 6% at +950 ms. Probabilities are computed over the labeled trials only. Ending the window at +800 ms instead gives the same reliabilities (0.98, 0.29, 0.66).
- **A few extreme subjects:** five subjects have overall recovery values of 0.45–0.58 (median 0.10). They raise the Pearson-based reliability; the rank-based value (0.81) still passes.
- **Error-effect reliability passes, but only marginally** (0.68; rank-based 0.57), so it is less dependable than the overall measure.

## Files

| File | Columns | Use |
|---|---|---|
| `results/features/recovery.csv` | `simon_recovery_correct`, `simon_recovery_err`, `simon_recovery_correct_matched`, `simon_recovery_error_effect`, counts `simon_recovery_n_correct`, `simon_recovery_n_err` | **Use** `simon_recovery_correct` and `simon_recovery_error_effect` for RQ3; the others are their components and bookkeeping |
| `results/features/recovery_conflict_unreliable.csv` | `simon_recovery_con`, `simon_recovery_inc`, `simon_recovery_conflict`, counts `simon_recovery_n_con`, `simon_recovery_n_inc` | **Do not use.** The conflict effect failed reliability (split-half r = 0.14, Spearman–Brown 0.24; rank-based 0.25; threshold 0.5). Kept only as a record of the pilot |

Both files have one row per `subject_id` and no clinical columns.
