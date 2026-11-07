"""B3: Rest-to-Simon microstate reorganization (K=5).

Run from the code folder: python B3.py
Inputs (data/): microstates_task.csv, subjects.csv
Output (results/): reorganization.csv  (one row per subject_id)

Definition of the reorganization index (reorganization_jsd):
    For each subject, take the Rest and Simon transition-probability
    distributions out of each source state s (the four off-diagonal
    probabilities P(s -> t), t != s; self-transitions excluded, rows sum
    to 1). Compute the Jensen-Shannon divergence (log base 2, range 0-1)
    between Rest and Simon for each source state, and average over the
    five source states A-E with equal weight:

        reorganization_jsd = (1/5) * sum_s JSD_2( P_rest(s -> .), P_simon(s -> .) )

Implementation: the 20 off-diagonal probabilities are concatenated and
each divided by 5, and one JSD (log 2) is taken. Because the five rows
have equal weight, this equals the mean of the five per-state JSDs
exactly (the KL terms decompose row by row).

Also saved: delta_cov_<state> = simon_cov - rest_cov (percentage points).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STATES = 'ABCDE'
# Row-major order: source A (to B,C,D,E), source B (to A,C,D,E), ...
EDGES = [(a, b) for a in STATES for b in STATES if a != b]


def transition_distribution(row, task):
    values = np.array([getattr(row, f'{task}_trans_{a}_to_{b}') for a, b in EDGES], dtype=float)
    matrix = values.reshape(5, 4)
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError(f'Invalid {task} transition probabilities for {row.subject_id}')
    if not np.allclose(matrix.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError(f'Transition rows not normalized for {row.subject_id} ({task})')
    return values / 5.0


def make_reorganization(metrics, subjects):
    # Exclude only explicit Simon trial label = exclude; keep ok and partial.
    eligibility = subjects[['subject_id', 'trials_simon']].copy()
    frame = metrics.merge(eligibility, on='subject_id', how='left', validate='one_to_one')
    required = ([f'{task}_cov_{s}' for task in ('rest', 'simon') for s in STATES]
                + [f'{task}_trans_{a}_to_{b}' for task in ('rest', 'simon') for a, b in EDGES])
    frame = frame.loc[frame.trials_simon.ne('exclude')].dropna(subset=required).copy()
    rows = []
    for row in frame.itertuples(index=False):
        rest = transition_distribution(row, 'rest')
        simon = transition_distribution(row, 'simon')
        # scipy returns sqrt(JSD); square to obtain divergence.
        reorg_jsd = float(jensenshannon(rest, simon, base=2.0) ** 2)
        result = {'subject_id': row.subject_id, 'reorganization_jsd': reorg_jsd}
        for state in STATES:
            result[f'delta_cov_{state}'] = getattr(row, f'simon_cov_{state}') - getattr(row, f'rest_cov_{state}')
        rows.append(result)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    microstates = pd.read_csv(DATA_DIR / "microstates_task.csv")
    subjects = pd.read_csv(DATA_DIR / "subjects.csv")

    result = make_reorganization(microstates, subjects)
    result.to_csv(RESULTS_DIR / "reorganization.csv", index=False)

    print('Subjects:', len(result))
    print(result.describe().round(4))
    print('Coverage delta row-sum maximum |error|:',
          result.filter(like='delta_cov_').sum(axis=1).abs().max())