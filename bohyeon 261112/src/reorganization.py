""" B3: Rest-to-Simon microstate reorganization (K states, detected from the CSV).

Run: python src/reorganization.py --micro microstate_A.csv [--tag A]
     (normally called by src/run_pipeline.py)
Inputs (data/): the run's microstate CSV, subjects.csv
Output: results/<tag>/reorganization_<tag>.csv  (one row per subject_id)

Definition of the reorganization index (reorganization_jsd):
    For each subject, take the Rest and Simon transition-probability
    distributions out of each source state s (the K-1 off-diagonal
    probabilities P(s -> t), t != s; self-transitions excluded, rows sum
    to 1). Compute the Jensen-Shannon divergence (log base 2, range 0-1)
    between Rest and Simon for each source state, and average over the
    K source states with equal weight:

        reorganization_jsd = (1/K) * sum_s JSD_2( P_rest(s -> .), P_simon(s -> .) )

Implementation: the K(K-1) off-diagonal probabilities are concatenated and
each divided by K, and one JSD (log 2) is taken. Because the K rows have
equal weight, this equals the mean of the K per-state JSDs exactly
(the KL terms decompose row by row). For K = 5 this is identical to the
original script.

Also saved: delta_cov_<state> = simon_cov - rest_cov (percentage points).
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from pipeline_io import (DATA_DIR, SUBJECTS_FILE, parse_args, load_microstates,
                         detect_states, save_run_config)


def make_edges(states):
    # Row-major order: source A (to B,C,...), source B (to A,C,...), ...
    return [(a, b) for a in states for b in states if a != b]


def transition_distribution(row, task, edges, k):
    values = np.array([getattr(row, f'{task}_trans_{a}_to_{b}') for a, b in edges], dtype=float)
    matrix = values.reshape(k, k - 1)
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError(f'Invalid {task} transition probabilities for {row.subject_id}')
    if not np.allclose(matrix.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError(f'Transition rows not normalized for {row.subject_id} ({task})')
    return values / k


def make_reorganization(metrics, subjects, states, simon_qc=True):
    edges = make_edges(states)
    k = len(states)
    eligibility = subjects[['subject_id', 'trials_simon']].copy()
    frame = metrics.merge(eligibility, on='subject_id', how='left', validate='one_to_one')
    required = ([f'{task}_cov_{s}' for task in ('rest', 'simon') for s in states]
                + [f'{task}_trans_{a}_to_{b}' for task in ('rest', 'simon') for a, b in edges])
    if simon_qc:
        # Exclude only explicit Simon trial label = exclude; keep ok and partial.
        frame = frame.loc[frame.trials_simon.ne('exclude')]
    frame = frame.dropna(subset=required).copy()
    rows = []
    for row in frame.itertuples(index=False):
        rest = transition_distribution(row, 'rest', edges, k)
        simon = transition_distribution(row, 'simon', edges, k)
        # scipy returns sqrt(JSD); square to obtain divergence.
        reorg_jsd = float(jensenshannon(rest, simon, base=2.0) ** 2)
        result = {'subject_id': row.subject_id, 'reorganization_jsd': reorg_jsd}
        for state in states:
            result[f'delta_cov_{state}'] = getattr(row, f'simon_cov_{state}') - getattr(row, f'rest_cov_{state}')
        rows.append(result)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    ctx = parse_args(__doc__)
    microstates = load_microstates(ctx)
    states = detect_states(microstates)
    subjects = pd.read_csv(DATA_DIR / SUBJECTS_FILE)

    result = make_reorganization(microstates, subjects, states, simon_qc=ctx.simon_qc)
    result.to_csv(ctx.csv("reorganization"), index=False)
    save_run_config(ctx, states, 'B3_reorganization')

    print('Subjects:', len(result))
    print(result.describe().round(4))
    print('Coverage delta row-sum maximum |error|:',
          result.filter(like='delta_cov_').sum(axis=1).abs().max())