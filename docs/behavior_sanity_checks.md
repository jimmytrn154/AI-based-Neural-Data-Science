# Simon behavior: sanity checks against the Dataset Guide

For the sanity-check slide (task M2). Computed Oct 23 by Chuong with `notebooks/chuong_behavior.py`; per-subject values are in `results/features/behavior.csv`.

**Sample:** 141 subjects with usable Simon trial labels (`trials_simon` = ok or partial): 94 PD, 47 controls. Each value is averaged within each subject first, then across subjects.

## Reference values

| Measure | Dataset Guide | Ours | Difference |
|---|---|---|---|
| RT, congruent (correct trials) | 514.7 ms | 514.4 ms | −0.3 ms |
| RT, incongruent (correct trials) | 557.4 ms | 557.4 ms | 0.0 ms |
| Error rate, congruent | 6.4% | 6.4% | 0.0 points |
| Error rate, incongruent | 14.8% | 14.8% | 0.0 points |

**All four reference values are reproduced.** The error rates match when calculated as errors divided by all trials, which is how the Dataset Guide calculates them.

`behavior.csv` uses a slightly different error rate: errors divided by trials that got a response, leaving out no-response trials. That version gives 6.8% (congruent) and 16.3% (incongruent). For the slide, use the Dataset Guide definition above.

## Simon effect by group (descriptive)

| | PD (n = 94) | Guide PD | Control (n = 47) | Guide control |
|---|---|---|---|---|
| Simon effect on RT (incongruent − congruent) | +47 ms | +50 ms | +35 ms | +34 ms |
| Error rate, congruent | 7.4% | 7.4% | 4.5% | 4.5% |
| Error rate, incongruent | 18.0% | 17.6% | 8.5% | 8.4% |

Error rates in this table use the Dataset Guide definition (errors / all trials). Both groups show the Simon effect. PD participants make about twice as many incongruent errors as controls, as in the Dataset Guide.

These are descriptive means only; the group comparisons with statistics are Bohyeon's (task B6).

## Suggested slide message

"Behavior reproduces the Dataset Guide: RTs within 0.3 ms, error rates exact; PD participants make about twice as many errors on conflict trials."
