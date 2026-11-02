r"""S5: paired cluster permutation test of Simon state probabilities.

Save alongside seonwoo_state_probabilities.py and seonwoo_erp_segmentation.py
in notebooks/. Run from the project root:
    .\.venv\Scripts\python.exe -X utf8 notebooks/seonwoo_conflict_stats.py

Optional software/subset check, using completed S2--S4 subset outputs:
    python -X utf8 notebooks/seonwoo_conflict_stats.py --subjects sub-P001 sub-C004

Uses saved S4 probabilities; no EEG filtering or back-fitting is repeated.
Contrast: incongruent minus congruent, one observation per participant.
Tests the configured analysis window, including its prestimulus portion.
Clusters connect adjacent times WITHIN a state, never different state labels.
The permutation maximum covers ALL states, times and both signs together.
Thus p_fwer is already familywise corrected, not an uncorrected pointwise p.
Missing probabilities cause a documented whole-participant exclusion for S5;
they are never filled with zeros. The S4 files remain unchanged.

References:
https://mne.tools/stable/generated/mne.stats.permutation_cluster_1samp_test.html
https://mne.tools/stable/generated/mne.stats.combine_adjacency.html
https://www.fieldtriptoolbox.org/faq/stats/clusterstats_interpretation/
"""

import argparse
import inspect
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src import io as pio
from seonwoo_erp_segmentation import choose_subjects
from seonwoo_state_probabilities import (
    digest_file, mean_sem, pipeline_signature, save_npz_atomic, subject_signature,
)

CLUSTER_COLUMNS = [
    "cluster_id", "state", "first_sample_ms", "last_sample_ms", "n_samples",
    "sample_bin_span_ms", "direction", "absolute_t_sum", "p_fwer",
    "passes_cluster_threshold", "mean_inc_minus_con_probability",
    "mean_inc_minus_con_percentage_points", "n_subjects",
]


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read_audit(path):
    table = pd.read_csv(path)
    check("subject_id" in table and "status" in table, f"Invalid audit: {path}")
    check(not table["subject_id"].isna().any(), f"Missing subject ID: {path}")
    check(not table["subject_id"].duplicated().any(), f"Duplicate subject ID: {path}")
    return table.set_index("subject_id")


def load_s4(cfg, requested, suffix):
    stem = f"seonwoo_state_probabilities{suffix}"
    path = pio.results_dir(cfg, "erp", stem + ".npz")
    if not path.is_file():
        raise FileNotFoundError(f"Complete S4 first. Missing: {path}")
    with np.load(path, allow_pickle=False) as file:
        data = {key: file[key] for key in file.files}
    metadata = json.loads(str(data["metadata_json"].item()))
    check(metadata["schema"] == "seonwoo_s4_v1", "Unknown S4 file schema.")
    expected_mode = "subset_check" if requested else "full_cohort"
    check(metadata["mode"] == expected_mode, "Full-cohort/subset file mismatch.")
    check(metadata["probability_denominator"] == "labeled_trials_at_each_timepoint",
          "Unexpected S4 probability denominator.")
    check(metadata["baseline_subtracted"] is False, "Unexpected probability baseline.")
    check(data["condition_codes"].tolist() == [1, 2], "Unexpected condition order.")
    check(data["condition_names"].tolist() == ["Congruent", "Incongruent"],
          "Unexpected condition names.")
    ids = data["subject_id"].tolist()
    states = data["state_names"].tolist()
    k = int(cfg["microstates"]["primary_k"])
    check(states == cfg["microstates"]["confirmed_labels"][k],
          "S4 state labels differ from the current config.")
    check(len(ids) == len(set(ids)), "Duplicate subjects in S4 arrays.")
    times = np.asarray(data["times_ms"], dtype=float)
    p = np.asarray(data["probabilities"], dtype=float)
    p_all = np.asarray(data["probabilities_all_trials"], dtype=float)
    nl = data["n_labeled"]
    nt = data["n_trials"]
    check(times.ndim == 1 and len(times) > 1 and np.isfinite(times).all(),
          "Invalid S4 time axis.")
    lo, hi = cfg["simon_conditions"]["analysis_window_ms"]
    dt = 1000.0 / cfg["preprocessing"]["resample_hz"]
    check(np.allclose(np.diff(times), dt) and times[0] >= lo - 1e-6
          and times[-1] <= hi + 1e-6 and times[0] < lo + dt
          and times[-1] > hi - dt, "S4 time grid/window differs from config.")
    check(metadata["analysis_window_ms"] == [lo, hi], "S4 window metadata mismatch.")
    check(p.shape == (len(ids), 2, k, len(times)), "S4 probability shape mismatch.")
    check(p_all.shape == p.shape and nl.shape == (len(ids), 2, len(times))
          and nt.shape == (len(ids), 2), "S4 denominator shape mismatch.")
    check(np.isfinite(nt).all() and np.isfinite(nl).all()
          and (nt == np.floor(nt)).all() and (nl == np.floor(nl)).all(),
          "Invalid S4 trial counts.")
    check((nt >= cfg["subjects"]["simon"]["min_trials_per_condition"]).all()
          and (nl >= 0).all() and (nl <= nt[:, :, None]).all(),
          "S4 trial counts violate the configured rules.")
    expected_valid = np.broadcast_to((nl > 0)[:, :, None, :], p.shape)
    check(np.array_equal(np.isfinite(p), expected_valid),
          "S4 missing probabilities do not match labeled-trial denominators.")
    check(np.isnan(p[~expected_valid]).all(), "Infinite S4 probability encountered.")
    check(((p[expected_valid] >= 0) & (p[expected_valid] <= 1)).all(),
          "S4 probability is outside [0, 1].")
    summed = np.where(expected_valid, p, 0).sum(axis=2)
    check(np.allclose(summed[nl > 0], 1), "S4 state probabilities do not sum to 1.")
    expected_all = np.where(expected_valid, p, 0) * (nl / nt[:, :, None])[:, :, None, :]
    check(np.isfinite(p_all).all() and np.allclose(p_all, expected_all),
          "S4 alternate probabilities do not agree with denominator counts.")

    # Detect stale aggregates, partial reruns, and changed cohort/config inputs.
    expected_ids, _ = choose_subjects(cfg, states, requested, suffix)
    check(ids == expected_ids, "S4 cohort/order differs from current S2 rules. Rerun S4.")
    s3 = read_audit(pio.results_dir(cfg, "features",
                                   f"seonwoo_erp_segmentation{suffix}_qc.csv"))
    finished = {"included", "excluded_subject_rule", "excluded_trial_status",
                "excluded_recording_qc", "excluded_min_trials", "excluded_microstate_qc"}
    check(s3["status"].isin(finished).all(), "S3 audit contains unfinished/error rows.")
    check(set(s3.index[s3["status"].eq("included")]) == set(ids),
          "S3/S4 included cohorts differ.")
    s4 = read_audit(pio.results_dir(cfg, "features", stem + "_qc.csv"))
    check(set(s4.index) == set(ids) and s4["status"].eq("included").all(),
          "S4 is incomplete: its audit does not cover every saved subject.")
    signature = pipeline_signature(cfg, k)
    check(metadata["pipeline_signature"] == signature,
          "Config/templates/S3-S4 code/packages changed since S4. Rerun S4 first.")
    for i, sid in enumerate(ids):
        counts = [int(s3.loc[sid, "n_congruent"]), int(s3.loc[sid, "n_incongruent"])]
        check(nt[i].tolist() == counts, f"{sid}: S3/S4 trial counts differ.")
        check([int(s4.loc[sid, "simon_n_con"]), int(s4.loc[sid, "simon_n_inc"])]
              == counts, f"{sid}: S4 audit counts differ.")
        cache = pio.results_dir(cfg, "labels", f"seonwoo_cue_K{k}", sid + ".npz")
        fingerprint = subject_signature(sid, cfg, signature, counts)
        with np.load(cache, allow_pickle=False) as saved:
            check(str(saved["fingerprint"].item()) == fingerprint,
                  f"{sid}: input files changed since S4. Rerun S4 first.")
            check(str(saved["subject_id"].item()) == sid
                  and np.array_equal(saved["times_ms"], times),
                  f"{sid}: S4 checkpoint identity/time mismatch.")
            for field in ("probabilities", "probabilities_all_trials", "n_labeled", "n_trials"):
                check(np.allclose(saved[field], data[field][i], equal_nan=True),
                      f"{sid}: aggregate differs from its S4 checkpoint. Finish rerunning S4.")
    return data, metadata, path


def paired_t(differences):
    """Ordinary paired t; all-zero columns have t=0 rather than undefined 0/0."""
    mean = differences.mean(axis=0)
    sem = differences.std(axis=0, ddof=1) / np.sqrt(len(differences))
    statistic = np.divide(mean, sem, out=np.zeros_like(mean), where=sem > 0)
    constant_nonzero = (sem == 0) & (mean != 0)
    statistic[constant_nonzero] = np.copysign(np.inf, mean[constant_nonzero])
    return statistic


def cluster_test(difference, n_permutations, cluster_forming_p, seed, verbose=True):
    check(difference.ndim == 3 and len(difference) >= 2,
          "Need at least 2 participants with complete paired curves.")
    check(np.isfinite(difference).all(), "Non-finite data entered the permutation test.")
    n, k, n_times = difference.shape
    threshold = float(student_t.ppf(1 - cluster_forming_p / 2, n - 1))
    x = difference.reshape(n, k * n_times).copy()
    observed = paired_t(x)
    if not np.any(np.abs(observed) > threshold):
        # No observed cluster means no candidate to test; MNE also skips in this case.
        return observed.reshape(k, n_times), [], np.empty(0), np.empty(0), threshold
    adjacency = mne.stats.combine_adjacency(np.zeros((k, k)), n_times)
    kwargs = dict(threshold=threshold, n_permutations=n_permutations, tail=0,
                  stat_fun=paired_t, adjacency=adjacency, n_jobs=1,
                  t_power=1, step_down_p=0, out_type="indices", check_disjoint=True,
                  buffer_size=None, verbose=verbose)
    # MNE 1.13 added rng; older project environments use seed. Record the API below.
    parameters = inspect.signature(mne.stats.permutation_cluster_1samp_test).parameters
    kwargs["rng" if "rng" in parameters else "seed"] = seed
    t_obs, clusters, p_values, h0 = mne.stats.permutation_cluster_1samp_test(x, **kwargs)
    return t_obs.reshape(k, n_times), clusters, p_values, h0, threshold


def cluster_table(clusters, p_values, t_obs, difference, times, states, alpha):
    k, n_times = t_obs.shape
    membership = np.zeros((k, n_times), dtype=np.int32)
    p_map = np.full((k, n_times), np.nan)
    rows = []
    dt = float(times[1] - times[0])
    for number, (cluster, p_value) in enumerate(zip(clusters, p_values), start=1):
        flat = np.asarray(cluster[0], dtype=int)
        state_indices, time_indices = np.unravel_index(flat, (k, n_times))
        check(len(np.unique(state_indices)) == 1, "A cluster crossed state labels.")
        s = int(state_indices[0])
        ti = np.sort(time_indices)
        check(np.all(np.diff(ti) == 1), "A cluster contains disconnected time samples.")
        sign = np.sign(t_obs[s, ti[0]])
        check(np.all(np.sign(t_obs[s, ti]) == sign), "A cluster crossed t-statistic signs.")
        membership[s, ti] = number
        p_map[s, ti] = p_value
        effect = float(difference[:, s, ti].mean())
        rows.append({
            "cluster_id": number, "state": states[s],
            "first_sample_ms": float(times[ti[0]]),
            "last_sample_ms": float(times[ti[-1]]), "n_samples": len(ti),
            "sample_bin_span_ms": len(ti) * dt,
            "direction": "inc>con" if sign > 0 else "inc<con",
            "absolute_t_sum": float(np.abs(t_obs[s, ti]).sum()),
            "p_fwer": float(p_value), "passes_cluster_threshold": bool(p_value < alpha),
            "mean_inc_minus_con_probability": effect,
            "mean_inc_minus_con_percentage_points": 100 * effect,
            "n_subjects": len(difference),
        })
    return pd.DataFrame(rows, columns=CLUSTER_COLUMNS), membership, p_map


def make_figure(times, probabilities, states, clusters, alpha, n_input, path, subset):
    mean, sem, _ = mean_sem(probabilities)
    dm, ds, _ = mean_sem(probabilities[:, 1] - probabilities[:, 0])
    k = len(states)
    fig, axes = plt.subplots(k, 2, figsize=(13.0, 10.6), sharex=True,
                             gridspec_kw={"width_ratios": [1.25, 1]}, squeeze=False)
    y_top = min(1, max(0.2, float(np.max(mean + sem)) * 1.12))
    d_top = max(1.0, float(np.max(np.abs(dm) + ds) * 100) * 1.25)
    dt = float(times[1] - times[0])
    for s, state in enumerate(states):
        left, right = axes[s]
        for c, (name, color) in enumerate((("Congruent", "#2463A6"),
                                            ("Incongruent", "#D65C3A"))):
            left.plot(times, mean[c, s], color=color, lw=1.4, label=name)
            left.fill_between(times, mean[c, s] - sem[c, s], mean[c, s] + sem[c, s],
                              color=color, alpha=0.16)
        left.set_ylim(0, y_top)
        left.set_ylabel(f"P({state})")
        right.plot(times, dm[s] * 100, color="#634A88", lw=1.5)
        right.fill_between(times, (dm[s] - ds[s]) * 100, (dm[s] + ds[s]) * 100,
                           color="#634A88", alpha=0.17)
        right.axhline(0, color="0.55", lw=0.8)
        right.set_ylim(-d_top, d_top)
        right.set_ylabel("Difference (pp)")
        supported = clusters.loc[(clusters["state"] == state)
                                 & clusters["passes_cluster_threshold"].eq(True)]
        for ax in (left, right):
            ax.axvline(0, color="0.55", ls="--", lw=0.8)
            ax.set_xlim(times[0], times[-1])
            ax.grid(axis="y", alpha=0.15)
            ax.spines[["top", "right"]].set_visible(False)
            for row in supported.itertuples():
                start = max(times[0], row.first_sample_ms - dt / 2)
                stop = min(times[-1], row.last_sample_ms + dt / 2)
                ax.plot([start, stop], [0.045, 0.045], transform=ax.get_xaxis_transform(),
                        lw=4, color="#202020", solid_capstyle="butt")
    axes[0, 0].set_title("Condition probabilities", fontsize=11)
    axes[0, 1].set_title("Incongruent minus congruent", fontsize=11)
    axes[0, 0].legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")
    for ax in axes[-1]:
        ax.set_xlabel("Time from cue onset (ms)")
    title = "Simon conflict | paired cluster permutation test"
    if subset:
        title = "SUBSET CHECK ONLY | " + title
    fig.suptitle(title + f" | N = {len(probabilities)}", fontsize=13, y=0.98)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.91, bottom=0.13,
                        hspace=0.22, wspace=0.25)
    n_supported = int(clusters["passes_cluster_threshold"].sum()) if len(clusters) else 0
    fig.text(0.5, 0.035,
             f"Means ±1 SEM; complete paired participants: {len(probabilities)}/{n_input}. "
             "pp = percentage points.\n"
             f"Black bars: membership in clusters with p_FWER < {alpha:g} "
             f"({n_supported} clusters); correction covers all states, tested times and both signs.\n"
             "Cluster bounds do not establish the exact onset or offset of an effect.",
             ha="center", fontsize=8.5, color="#444444")
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", help="Subset check; requires S2--S4 subset files.")
    args = parser.parse_args()
    suffix = "_subset" if args.subjects else ""
    cfg = pio.load_config()
    settings = cfg["simon_conditions"]["permutation_test"]
    check(settings["type"] == "paired cluster-based", "Unexpected configured test type.")
    n_permutations = int(settings["n_permutations"])
    forming_p = float(settings["cluster_forming_p"])
    alpha = float(cfg["statistics"]["alpha"])
    seed = int(cfg["random_seed"])
    check(n_permutations >= 2 and 0 < forming_p < 1 and 0 < alpha < 1,
          "Invalid permutation-test settings.")
    print("S5: checking saved S4 probabilities, cohort and checkpoints...", flush=True)
    data, s4_metadata, input_path = load_s4(cfg, args.subjects, suffix)
    ids, times, states = data["subject_id"], data["times_ms"], data["state_names"].tolist()
    p = data["probabilities"]
    complete = np.isfinite(p).all(axis=(1, 2, 3))
    stats_dir = pio.results_dir(cfg, "stats")
    figure_dir = pio.results_dir(cfg, "figures")
    stats_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    stem = f"seonwoo_conflict{suffix}"
    subject_audit = pd.DataFrame({
        "subject_id": ids,
        "status": np.where(complete, "included", "excluded_nonfinite_probability"),
        "missing_probability_values": (~np.isfinite(p)).sum(axis=(1, 2, 3)),
        "simon_n_con": data["n_trials"][:, 0], "simon_n_inc": data["n_trials"][:, 1],
    })
    subject_audit.to_csv(stats_dir / f"{stem}_subjects.csv", index=False)
    n = int(complete.sum())
    check(n >= 2, "Fewer than 2 complete participants; see the S5 subject audit.")
    paired = p[complete]
    difference = paired[:, 1] - paired[:, 0]
    print(f"Complete paired subjects: {n}/{len(ids)}; excluded for missing data: {len(ids)-n}")
    print(f"Contrast: incongruent - congruent; states {states}; "
          f"window {times[0]:g} to {times[-1]:g} ms.")
    print(f"Requested permutations: {n_permutations}; two-sided; "
          "correction across ALL states and times.", flush=True)
    t_obs, clusters, p_values, h0, threshold = cluster_test(
        difference, n_permutations, forming_p, seed
    )
    table, membership, p_map = cluster_table(
        clusters, p_values, t_obs, difference, times, states, alpha
    )
    table_path = stats_dir / f"{stem}_clusters.csv"
    table.to_csv(table_path, index=False)
    mean, sem, _ = mean_sem(paired)
    dm, ds, _ = mean_sem(difference)
    curves = []
    for s, state in enumerate(states):
        curves.append(pd.DataFrame({
            "state": state, "time_ms": times, "n_paired_subjects": n,
            "con_mean": mean[0, s], "inc_mean": mean[1, s],
            "con_sem": sem[0, s], "inc_sem": sem[1, s],
            "inc_minus_con_mean": dm[s], "inc_minus_con_sem": ds[s],
            "paired_t": t_obs[s], "cluster_id": membership[s],
            "cluster_p_fwer": p_map[s],
            "member_of_supported_cluster": p_map[s] < alpha,
        }))
    pd.concat(curves, ignore_index=True).to_csv(stats_dir / f"{stem}_timecourse.csv", index=False)
    n_supported = int(table["passes_cluster_threshold"].sum()) if len(table) else 0
    exact = n_permutations >= 2 ** (n - 1)
    metadata = {
        "schema": "seonwoo_s5_v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": s4_metadata["mode"], "input_s4_sha256": digest_file(input_path),
        "input_s4_pipeline_signature": s4_metadata["pipeline_signature"],
        "script_sha256": digest_file(__file__), "mne_version": mne.__version__,
        "contrast": "incongruent_minus_congruent", "probability_denominator":
        s4_metadata["probability_denominator"], "n_s4": len(ids), "n_paired": n,
        "missing_rule": "complete participant across both conditions/all states/entire window",
        "window_ms": [float(times[0]), float(times[-1])],
        "tail": 0, "cluster_forming_p": forming_p, "cluster_forming_t": threshold,
        "alpha": alpha, "n_permutations_requested": n_permutations,
        "n_null_statistics": len(h0), "exact_if_clusters_present": exact, "seed": seed,
        "random_api": "rng" if "rng" in inspect.signature(
            mne.stats.permutation_cluster_1samp_test).parameters else "seed",
        "adjacency": "adjacent times within each state; no cross-state connections",
        "correction": "maximum absolute cluster t-sum across all states/times/both signs",
        "n_clusters": len(table), "n_clusters_passing_threshold": n_supported,
        "team_general_fdr_setting": cfg["statistics"].get("fdr"),
    }
    save_npz_atomic(stats_dir / f"{stem}_test.npz", subject_id=ids[complete],
                    times_ms=times, state_names=np.asarray(states), paired_t=t_obs,
                    cluster_id=membership, cluster_p_fwer=p_map, null_max_cluster=h0,
                    metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)))
    figure_path = figure_dir / f"seonwoo_conflict_timecourse{suffix}.png"
    make_figure(times, paired, states, table, alpha, len(ids), figure_path, bool(args.subjects))
    lines = [
        "S5: Simon congruent versus incongruent state probabilities",
        "SUBSET CHECK ONLY" if args.subjects else "Full eligible cohort",
        f"S4 subjects: {len(ids)}; S5 complete paired subjects: {n}; missing-data exclusions: {len(ids)-n}.",
        f"Window: {times[0]:g} to {times[-1]:g} ms, inclusive; states: {states}.",
        "Contrast: incongruent minus congruent. Each participant has equal weight.",
        "Primary probabilities condition on labeled trials; no probability baseline subtraction.",
        "Missing rule: exclude a participant from S5 if any required probability is non-finite.",
        "Paired sign flips act on an entire participant's difference array, preserving state/time dependence.",
        "Assumption: independent participants and sign-exchangeable paired differences under the null.",
        f"Two-sided paired t; cluster-forming p={forming_p:g}, t threshold={threshold:.6g}, df={n-1}.",
        "All-zero difference columns have t=0; constant nonzero columns have signed infinite t.",
        "Clusters connect neighboring time samples within a state, with separate positive/negative clusters.",
        "Cluster statistic: absolute sum of t; null maximum over ALL states, times and both signs.",
        f"Requested permutations: {n_permutations}; actual null statistics: {len(h0)}; seed={seed}.",
        f"MNE random API: {metadata['random_api']} (random streams may differ across MNE versions).",
        "Exact enumeration replaces sampling when the number of distinct two-sided sign patterns is smaller.",
        "If no observed cluster exceeds the forming threshold, permutations are skipped and the table is empty.",
        f"Observed clusters: {len(table)}; clusters with familywise-corrected p < {alpha:g}: {n_supported}.",
        "S5 correction is max-cluster FWER across this state/time family. Do not label these p-values as FDR q-values.",
        "No second BH adjustment is applied to selected clusters; planned scalar-feature families still require the team's BH rule.",
        f"General config FDR note: {cfg['statistics'].get('fdr', 'not specified')}",
        "Cluster membership supports the condition-level test, not significance of individual time samples or precise onset/offset.",
        "Cluster-averaged differences in the table are descriptive, selected using the same data, not independent effect estimates.",
        "p >= alpha does not establish equivalence or absence of a condition difference.",
        "Figures show means +/-1 SEM over the SAME complete participants used in the test.",
        "S4 trial order, smoothing, cohort and template-QC assumptions are inherited unchanged.",
        "No PD/control or clinical association test is performed here; S6 conflict.csv is a separate next step.",
        f"Input SHA256: {metadata['input_s4_sha256']}",
        f"MNE version: {mne.__version__}",
        "Reference: https://mne.tools/stable/generated/mne.stats.permutation_cluster_1samp_test.html",
        "Reference: https://mne.tools/stable/generated/mne.stats.combine_adjacency.html",
        "Interpretation: https://www.fieldtriptoolbox.org/faq/stats/clusterstats_interpretation/",
    ]
    summary_path = stats_dir / f"{stem}_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nClusters passing corrected p < {alpha:g}: {n_supported}")
    if not len(table):
        print("No observed clusters exceeded the forming threshold; no permutations needed.")
    elif len(h0) < n_permutations:
        print(f"Exact small-sample test: {len(h0)} distinct null statistics.")
    if args.subjects:
        print("SUBSET CHECK: these outputs are not the full-cohort project result.")
    print(f"Saved figure: {figure_path}\nSaved cluster table: {table_path}\nSaved summary: {summary_path}")
    print("DONE: conflict cluster test completed." if not args.subjects
          else "DONE: conflict cluster subset check completed.")


if __name__ == "__main__":
    main()
