# Slide content: Chuong, final presentation (Phase 2)

Draft for the shared final deck. One section per slide. Every number comes from the committed results; the source is in brackets after each figure. State E is called **task-enhanced**, never "task-specific".

---

## Slide 1. Is E-up just C-down redistributed? (P2.2)

- Coverage percentages sum to 100, so when C loses time, other states must gain. Question: does E gain **more** than its fair share?
- If C's lost time were shared out in proportion to Rest coverage, E would rise by about 0.6 points. It rises by 1.8, so about two-thirds is excess.
- E's share of the non-C time rises: +1.83 points [+0.35, +3.35], FDR q 0.0035 (N = 139).
- In log-ratio terms (which remove the sum-to-100 constraint), E shows the largest increase of the five states (CLR +0.18, q 0.0003) and rises relative to B and D.
- Without the competition between states: the average fit of every sample to the E map rises (mean |r| +0.007, q 0.012), while the fit to C, B and D falls.
- E samples in Simon are not more C-like (C as runner-up 28.8% → 28.1%, q 0.035), so E-up is not mainly C frames relabelled as E.
- Same conclusions with blink-cleaned Rest. Caveat: the effects are small, and the E fit does not rise significantly at GFP peaks.

**Figure:** `results/figures/chuong_p2_2_probability_dependence.png` and `results/figures/chuong_p2_2_part2.png` [`PROGRESS.md`, P2.2 notes].

**Speaking note:** A drop in one state forces a rise somewhere else, so we tested E against that constraint, with the rules fixed before computing. E-up is not explained by C-down alone. Looking at each sample's fit to every map, E-like activity increases in general rather than C frames tipping over into E.

---

## Slide 2. Are the results driven by blinks? (P2.7)

- 33 of the 139 Rest recordings had no eye component removed by ICA. We cleaned them (0.5-s windows with Fp1/Fp2 > 100 µV excluded; rule fixed in advance) and reran everything that uses Rest.
- **Rest → Simon change:** E holds (+1.67 points, p 0.0017); C weakens (−2.26, p 0.072). A holds.
- **Group results using Rest:** 39 of 40 RQ1 mixed-model terms and 106 of 108 RQ3 tests keep their conclusion. PD still has less E at rest (−4.07 points, q 0.010).
- **Conflict contrasts without the 9 Simon recordings with no eye component removed:** P2-B7 6 of 6 hold. P2-B8: the main result holds (N2 Δ(C) PD − control +2.47 points, q 0.010); three borderline results (q 0.041–0.044) move to q 0.059–0.069 with unchanged estimates (fewer subjects, not blinks).

**Figure:** none yet; table from `results/features/p2_7_blink_sensitivity.csv` [`PROGRESS.md`, P2.7 notes].

**Speaking note:** The E increase and the group results do not depend on blink artefacts. The C decrease is partly sensitive, so we present it with that caveat.

---

## Slide 3. Recovery after the response: go/no-go (P2.9)

- Pilot question: after the keypress, how quickly do the brain-state probabilities return to each person's pre-stimulus baseline? Measure: distance from baseline, averaged over 200–1000 ms after the keypress (140 subjects).
- Rule fixed in advance: keep a measure only if split-half reliability (Spearman–Brown) ≥ 0.5.
  - Overall recovery: 0.97 → **go**.
  - Errors vs. size-matched correct trials: 0.68 → **go, marginal** (after errors the state distribution stays further from baseline; median +0.020, higher in 67% of subjects).
  - Incongruent vs. congruent: 0.24 → **no-go**, not interpreted.
- Not yet tested against group or cognition; exploratory only.

**Figure:** `results/figures/chuong_recovery.png` [`docs/recovery_pilot.md`].

**Speaking note:** We set a reliability bar before looking. Two measures pass and can be used in later analyses as exploratory features; the conflict version fails, so we drop it rather than interpret it.
