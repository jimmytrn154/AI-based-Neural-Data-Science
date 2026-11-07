"""B6: RQ3 clinical associations, including conflict effect (N2/theta).

Run from the code folder: python B6.py
Inputs: results/reorganization.csv (B3); data/conflict.csv,
data/microstates_task.csv, data/behavior.csv, data/subjects.csv;
B1 functions from data/b1_outputs/stats.py.

Eighteen subject-level EEG features (config primary features):
reorganization_jsd, delta_cov_A-E, conflict_n2_l2, conflict_theta_l2,
rest_meandur_A-E, simon_meandur_A-E. Clinical variables are never saved to
the EEG feature table.
Conflict magnitudes are the Euclidean norm of signed state-wise conflict
probability differences, per window.
FDR: Benjamini-Hochberg over ONE RQ3 family (all group tests and all
partial correlations, config statistics.fdr). Mann-Whitney U (PD vs
Control, unadjusted) is reported as a check next to the linear model.
EEG features are merged with an outer join, so a subject keeps every
feature it has: e.g. sub-P001 (Rest excluded, Simon valid) enters the
conflict and Simon duration analyses.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

import sys
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(DATA_DIR / "b1_outputs"))
from stats import fdr_bh, group_linear_model, partial_spearman

EEG_BASE = ['reorganization_jsd'] + [f'delta_cov_{s}' for s in 'ABCDE']
EEG_CONFLICT = ['conflict_n2_l2', 'conflict_theta_l2']
EEG_DURATION = ([f'rest_meandur_{s}' for s in 'ABCDE']
                + [f'simon_meandur_{s}' for s in 'ABCDE'])
EEG_FEATURES = EEG_BASE + EEG_CONFLICT + EEG_DURATION
FAMILY = 'RQ3'   # one family for all RQ3 tests (config)
N_BOOT = 2000
SEED = 42


def locate(filename, folders=None):
    paths = [DATA_DIR / filename, RESULTS_DIR / filename]

    for path in paths:
        if path.is_file():
            return path

    raise FileNotFoundError(f"{filename}: searched {[str(p) for p in paths]}")

def load_unique(filename, folders):
    path = locate(filename, folders)
    table = pd.read_csv(path)
    if table['subject_id'].duplicated().any():
        raise ValueError(f'Duplicate subject_id in {path}')
    print(f'Loaded {filename}: {len(table)} rows ({path})')
    return table


def main():
    reorg = load_unique('reorganization.csv', ['b3_outputs', '', 'b6_outputs'])
    conflict = load_unique('conflict.csv', ['b6_outputs', '', 'b2_outputs'])
    behavior = load_unique('behavior.csv', ['', 'data_outputs', 'b6_outputs'])
    subjects = load_unique('subjects.csv', ['', 'data_outputs', 'b6_outputs'])
    micro = load_unique('microstates_task.csv', ['', 'data_outputs'])
    # Reference table is a summary of Wilcoxon tests, NOT subject-level input.
    try:
        rq1 = pd.read_csv(locate('rq1_stats.csv', ['b2_outputs', '']))
        print(f'B2 reference loaded: {len(rq1)} tests; not merged as patient features')
    except FileNotFoundError:
        print('Optional B2 reference not found; B6 proceeds normally')

    # For each subject, five state-specific signed conflict differences form a vector.
    # Its L2 norm is nonnegative and captures the magnitude of a condition difference.
    for window in ['n2', 'theta']:
        columns = [f'simon_conflict_{window}_{s}' for s in 'ABCDE']
        conflict[f'conflict_{window}_l2'] = np.sqrt(
            np.square(conflict[columns].to_numpy(dtype=float)).sum(axis=1)
        )

    # Keep a one-subject-one-row EEG feature table, without clinical columns.
    # Outer joins: Rest-based features are missing for sub-P001, but its
    # conflict and Simon duration features are kept.
    eeg = reorg[['subject_id'] + EEG_BASE].merge(
        conflict[['subject_id'] + EEG_CONFLICT], on='subject_id',
        how='outer', validate='one_to_one'
    ).merge(
        micro[['subject_id'] + EEG_DURATION], on='subject_id',
        how='left', validate='one_to_one'
    )
    eeg.to_csv(RESULTS_DIR / 'b6_eeg_features.csv', index=False)

    data = eeg.merge(
        subjects[['subject_id', 'group', 'age', 'gender', 'moca', 'updrs', 'trials_simon']],
        on='subject_id', how='inner', validate='one_to_one'
    ).merge(
        behavior[['subject_id', 'simon_err_inc']],
        on='subject_id', how='left', validate='one_to_one'
    )
    data = data[data['trials_simon'].isin(['ok', 'partial'])].copy()
    print('Eligible EEG subjects:',len(data),'PD/control:',data.group.value_counts().to_dict())
    print('Available reorganization (Rest+Simon):', data['reorganization_jsd'].notna().sum())
    print('Available N2/theta conflict:',data[EEG_CONFLICT].notna().all(axis=1).sum(),
          '| includes sub-P001:', bool(data.loc[data.subject_id.eq('sub-P001'), EEG_CONFLICT].notna().all(axis=1).any()))
    print('Available Simon error:', data.simon_err_inc.notna().sum())

    group_rows = []
    for feature in EEG_FEATURES:
        output = group_linear_model(data, feature)
        # Mann-Whitney check (unadjusted, no covariates; not in the FDR family)
        pd_values = pd.to_numeric(data.loc[data.group.eq('PD'), feature], errors='coerce').dropna()
        hc_values = pd.to_numeric(data.loc[data.group.eq('Control'), feature], errors='coerce').dropna()
        mw = mannwhitneyu(pd_values, hc_values, alternative='two-sided')
        group_rows.append({'analysis':'PD_vs_Control', 'feature':feature,
                           'family':FAMILY, **output,
                           'mw_u': float(mw.statistic), 'mw_p': float(mw.pvalue),
                           'mw_n_pd': len(pd_values), 'mw_n_control': len(hc_values)})

    corr_rows = []
    for subgroup in ['all', 'PD']:
        subset = data if subgroup=='all' else data[data.group.eq('PD')]
        for outcome in ['moca', 'simon_err_inc']:
            for feature in EEG_FEATURES:
                output = partial_spearman(subset, feature, outcome,
                                          n_boot=N_BOOT, seed=SEED)
                corr_rows.append({'analysis':'partial_spearman', 'feature':feature,
                                  'outcome':outcome, 'subgroup':subgroup,
                                  'family':FAMILY, **output})
    subset = data[data.group.eq('PD')]
    for feature in EEG_FEATURES:
        output = partial_spearman(subset, feature, 'updrs',
                                  n_boot=N_BOOT, seed=SEED)
        corr_rows.append({'analysis':'partial_spearman', 'feature':feature,
                          'outcome':'updrs', 'subgroup':'PD',
                          'family':FAMILY, **output})

    # One BH correction over all RQ3 tests, then split back into the two tables
    rq3 = fdr_bh(pd.concat([pd.DataFrame(group_rows), pd.DataFrame(corr_rows)],
                           ignore_index=True))
    gr = rq3[rq3.analysis.eq('PD_vs_Control')].dropna(axis=1, how='all')
    cr = rq3[rq3.analysis.eq('partial_spearman')].dropna(axis=1, how='all')
    gr.to_csv(RESULTS_DIR / 'b6_group_stats.csv', index=False)
    cr.to_csv(RESULTS_DIR / 'b6_correlations.csv', index=False)
    rq3.to_csv(RESULTS_DIR / 'rq3_stats.csv', index=False)
    print('Number of tests:',len(gr)+len(cr),'FDR-significant:',int(gr.reject_fdr.sum()+cr.reject_fdr.sum()))
    print('\nGroup results:\n', gr[['feature','n','estimate_pd_minus_control','p_value','p_fdr','mw_p']].to_string(index=False))
    print('\nConflict vs cognition and behavior:\n', cr[cr.feature.isin(EEG_CONFLICT)][['feature','outcome','subgroup','n','rho','p_value','p_fdr']].to_string(index=False))

    # Correlation summary plot: compact heat map of rho, NOT significance.
    keys = [('moca','all'),('moca','PD'),('simon_err_inc','all'),('simon_err_inc','PD'),('updrs','PD')]
    heat = np.full((len(EEG_FEATURES),len(keys)),np.nan)
    for j,(outcome,subgroup) in enumerate(keys):
        for i,feature in enumerate(EEG_FEATURES):
            row = cr[(cr.feature==feature)&(cr.outcome==outcome)&(cr.subgroup==subgroup)]
            if len(row): heat[i,j]=row.iloc[0]['rho']
    fig,ax=plt.subplots(figsize=(10.8,8.5))
    im=ax.imshow(heat,cmap='RdBu_r',vmin=-0.4,vmax=0.4,aspect='auto')
    ax.set_xticks(np.arange(len(keys)), ['MoCA\nAll','MoCA\nPD','Simon error\nAll','Simon error\nPD','UPDRS\nPD'])
    ax.set_yticks(np.arange(len(EEG_FEATURES)),EEG_FEATURES)
    for i in range(len(EEG_FEATURES)):
        for j in range(len(keys)):
            if np.isfinite(heat[i,j]): ax.text(j,i,f'{heat[i,j]:+.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,label='Age-adjusted partial Spearman rho')
    ax.set_title('B6: EEG microstate features and clinical/behavioral scores')
    fig.tight_layout()
    fig.savefig(FIGURES_DIR/'bohyeon_b6_correlations.png',dpi=220)
    plt.close(fig)
    print('\nOutputs saved in', RESULTS_DIR, 'and', FIGURES_DIR)


if __name__=='__main__':
    main()