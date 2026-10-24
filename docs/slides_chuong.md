# Slide content: Chuong (task C6), progress presentation Oct 26

Draft for the shared deck. One section per slide. Every number comes from the committed results; the source is in brackets after each figure. State E is called **task-enhanced**, never "task-specific".

---

## Slide 1. Pipeline: from cleaned EEG to brain-state sequences

- 149 Rest and 141 Simon recordings after the pre-specified exclusions; band-pass 2–20 Hz, resampled to 250 Hz.
- Templates: global field power (GFP) peaks, the same number from Rest and Simon per subject (up to 3,000 each); polarity-invariant modified k-means, first per subject, then across subjects.
- Data-quality rule fixed in advance: a recording is excluded if the templates fit it unusually badly (robust z < −3). At K = 5 this removes 2 recordings (one Rest, one Simon).
- Back-fitting with the temporal smoothing documented by pycrostates; our faster implementation gives labels identical to pycrostates (0 mismatching samples, tested).

**Figure:** none yet. A simple flow diagram (EEG → GFP peaks → two-level clustering → QC → back-fitting → metrics) to draw.

**Speaking note:** All rules, including the exclusions, filters and the quality-control threshold, were written into a shared config file before we looked at any results. Templates come from both tasks equally, so the same states can be compared between rest and the Simon task.

---

## Slide 2. Choosing the number of states

- We fitted K = 4 to 7. Template GEV: 54.8%, 56.8%, 58.4%, 59.3% (gains +2.0, +1.6, +0.9 points).
- Our pre-specified rule averaged four cluster scores and **pointed to K = 4**. Two of the four (Davies–Bouldin, Dunn) turned out degenerate at every K with pycrostates' distance measure, so the rule did not work as intended. The two that worked (silhouette, Calinski–Harabasz) also prefer K = 4.
- At K = 4 the canonical fronto-central map D disappears: it merges into C, and the fourth map is a different, central map. **Only K = 5 recovers canonical A–D** (r ≥ 0.97 with our draft maps), with the central map as a separate state E.
- We chose K = 5 before any group or clinical analysis, and keep K = 4 as a sensitivity analysis.

**Figure:** `results/figures/chuong_templates_K4.png` above `results/figures/chuong_templates_K5.png` [decision log, Oct 23; `results/templates/k_selection.csv`].

**Speaking note:** We want to be transparent here: our own selection rule pointed to four states, and we overrode it. The reason is that half of the rule's criteria failed numerically, and at four states one of the standard maps simply isn't there; five states keep all four canonical maps plus one new one, and the gain in explained variance flattens after five.

---

## Slide 3. The five brain states (K = 5)

- Group templates from 139 subjects; they explain 56.8% of the variance at GFP peaks.
- A, B, C, D match the canonical microstate classes; **E** is a central maximum with opposite-sign temporal edges.
- Labels were assigned by correlation with our approved draft maps, then confirmed visually.
- Typical dynamics: states last about 100 ms (median segment 96 ms; mean duration 102–153 ms by state); coverage sums to 100% in every recording.

**Figure:** `results/figures/chuong_templates_K5.png` [`results/features/microstates_task.csv`].

**Speaking note:** These five maps are what every later analysis uses, in both tasks and in every condition. The durations of around 100 milliseconds are in the range usually reported for EEG microstates.

---

## Slide 4. Do the same states describe rest and the Simon task?

- We refitted the K = 5 templates from Rest data only and from Simon data only.
- All five maps reappear in both: |r| with the shared templates ≥ 0.98 for every state, in both tasks.
- So rest and task can be compared on one common set of states (RQ1, first part).
- Sensitivity: at K = 4 the Rest-only fit did not contain the central map, another reason to use K = 5.

**Figure:** `results/figures/chuong_templates_K5_rest.png` and `results/figures/chuong_templates_K5_simon.png` side by side [decision log; `docs/TEAM_TASKS.md`].

**Speaking note:** Before comparing rest with the task, we checked that the states themselves are the same. They are: the task changes how often and how long the states occur, not which states exist, which is what makes the rest-versus-task comparison meaningful.

---

## Slide 5. State E is task-enhanced

- E is present at rest and during the task: median coverage 11.0% at rest, 12.1% in Simon; median duration about 100 ms in both.
- In Simon trials (20 subjects, correct trials), P(E) rises from 0.11 before the stimulus to a peak of 0.16 at 452 ms, close to the central positivity in the ERP (~440 ms at Cz), and peaks again at the keypress (+16 ms), then drops.
- E is not the muscle-like map from our first, unfiltered draft (|r| = 0.09).
- Descriptive only so far; the formal rest-versus-Simon tests are in the RQ1 results.

**Figure:** `results/figures/chuong_E_timecourse.png` [`notebooks/chuong_state_E.py`].

**Speaking note:** E is not something that only appears during the task: it is there at rest, and it becomes more frequent around the response. That is why we call it task-enhanced rather than task-specific. Because the new map has strong temporal edges, we checked that it is not a muscle artifact, and it does not resemble the muscle map we saw in the first draft.

---

## Slide 6. Exploratory: recovery after the response

- Measure: how far the brain-state probabilities are from each person's pre-stimulus baseline, averaged over 200–1000 ms after the keypress.
- Kept only if reliable (split-half, Spearman–Brown ≥ 0.5, rule fixed in advance): overall recovery 0.97 (kept); errors vs. correct 0.68 (kept, marginal); incongruent vs. congruent 0.24 (**dropped**).
- After errors, the state distribution stays further from baseline than after correct trials of the same number (median difference +0.020; higher in 67% of subjects).
- Errors involve fewer trials, which alone inflates the distance, so they are compared with size-matched sets of correct trials.

**Figure:** `results/figures/chuong_recovery.png` [`docs/recovery_pilot.md`].

**Speaking note:** This was a pilot, and we set a reliability bar before running it. The overall recovery measure and the error effect pass; the conflict effect does not, so we drop it rather than interpret it. The error effect will go into the group and cognition analyses only as an exploratory measure.

---

## Slide 7. Next phases

- **Phase 2 (after today's feedback):** robustness of the main results to K = 4, no smoothing, the QC-excluded recordings and possible blink contamination (files already prepared); alpha control with spectral parameterization; PCA trajectories for the recovery pilot.
- **Phase 3:** does EEG add value beyond demographics and behavior for detecting cognitive impairment (MoCA < 26)? Calibrated, subject-level prediction.
- **Phase 4:** the same states and measures in the Oddball and Interval-timing tasks; comparison with a deep model (EEGNet).
- Open question for the professor: how to interpret state E (see our questions slide).

**Figure:** none (text slide).

**Speaking note:** Phase 1 gave us a validated pipeline, a common set of five states and reliable task-level measures. The sensitivity files for the robustness checks are already prepared, so Phase 2 can start as soon as we have your feedback.
