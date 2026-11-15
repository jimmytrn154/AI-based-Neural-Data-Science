"""B6: RQ3 clinical associations, including conflict effect (N2/theta).

Run: python src/B6_clinical_associations.py --micro microstate_A.csv --conflict conflict_A.csv [--tag A]
     (normally called by src/run_pipeline.py; reorganization.py, and
     optionally B2, must run first with the same tag)
Inputs: results/<tag>/reorganization_<tag>.csv (B3); data/<conflict CSV>,
data/<microstate CSV>, data/behavior.csv, data/subjects.csv;
B1 functions from src/stats.py.
Outputs: results/<tag>/b6_eeg_features_<tag>.csv, b6_group_stats_<tag>.csv,
b6_correlations_<tag>.csv, rq3_stats_<tag>.csv;
figures/<tag>/bohyeon_b6_correlations_<tag>.png

Subject-level EEG features (K microstates detected from the CSV):
reorganization_jsd, delta_cov_<S>, conflict_n2_l2, conflict_theta_l2,
rest_meandur_<S>, simon_meandur_<S>  (18 features for K = 5, 15 for K = 4).
Clinical variables are never saved to the EEG feature table.
Conflict magnitudes are the Euclidean norm of signed state-wise conflict
probability differences, per window (states taken from the conflict file).
FDR: Benjamini-Hochberg over ONE RQ3 family (all group tests and all
partial correlations, config statistics.fdr). Mann-Whitney U (PD vs
Control, unadjusted) is reported as a check next to the linear model.
EEG features are merged with an outer join, so a subject keeps every
feature it has: e.g. sub-P001 (Rest excluded, Simon valid) enters the
conflict and Simon duration analyses.
"""
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"

sys.path.insert(0, str(HERE))
from scipy.stats import mannwhitneyu
from stats import fdr_bh, group_linear_model, partial_spearman
from pipeline_io import (SUBJECTS_FILE, BEHAVIOR_FILE, parse_args, load_microstates,
                         detect_states, apply_simon_qc, save_run_config)

FAMILY = 'RQ3'   # one family for all RQ3 tests (config)
N_BOOT = 2000
SEED = 42
EEG_CONFLICT = ['conflict_n2_l2', 'conflict_theta_l2']


def locate(filename, folder):
    path = folder / filename
    if not path.is_file():
        raise FileNotFoundError(f"{filename}: not found in {folder}")
    return path


def load_unique(path):
    table = pd.read_csv(path)
    if table['subject_id'].duplicated().any():
        raise ValueError(f'Duplicate subject_id in {path}')
    print(f'Loaded {path.name}: {len(table)} rows ({path})')
    return table


def main():
    ctx = parse_args(__doc__)
    if not ctx.conflict_file:
        raise ValueError("B6 needs --conflict (e.g. conflict_A.csv)")
    RESULTS_DIR, FIGURES_DIR = ctx.results_dir, ctx.figures_dir

    # reorganization_<tag>.csv must come from THIS run's results folder
    reorg_path = ctx.csv('reorganization')
    if not reorg_path.is_file():
        raise FileNotFoundError(
            f'{reorg_path} not found: run src/reorganization.py --micro {ctx.micro_file} --tag {ctx.name} first')
    reorg = load_unique(reorg_path)
    conflict = load_unique(locate(ctx.conflict_file, DATA_DIR))
    behavior = load_unique(locate(BEHAVIOR_FILE, DATA_DIR))
    subjects = load_unique(locate(SUBJECTS_FILE, DATA_DIR))
    micro = load_microstates(ctx)
    states = detect_states(micro)

    EEG_BASE = ['reorganization_jsd'] + [f'delta_cov_{s}' for s in states]
    EEG_DURATION = ([f'rest_meandur_{s}' for s in states]
                    + [f'simon_meandur_{s}' for s in states])
    EEG_FEATURES = EEG_BASE + EEG_CONFLICT + EEG_DURATION
    missing = [c for c in EEG_BASE if c not in reorg.columns]
    if missing:
        raise ValueError(f'{reorg_path} lacks {missing}: re-run reorganization.py --micro {ctx.micro_file} --tag {ctx.name}')
    print(f'EEG features: {len(EEG_FEATURES)}')

    # Reference table is a summary of Wilcoxon tests, NOT subject-level input.
    rq1_path = ctx.csv('rq1_stats')
    if rq1_path.is_file():
        rq1 = pd.read_csv(rq1_path)
        print(f'B2 reference loaded: {len(rq1)} tests; not merged as patient features')
    else:
        print('Optional B2 reference not found; B6 proceeds normally')

    # For each subject, the state-specific signed conflict differences form a vector.
    # Its L2 norm is nonnegative and captures the magnitude of a condition difference.
    # States are those in the conflict file.
    for window in ['n2', 'theta']:
        columns = sorted(c for c in conflict.columns
                         if re.fullmatch(rf'simon_conflict_{window}_[A-Z]', c))
        if not columns:
            raise ValueError(f'No simon_conflict_{window}_<S> columns in {ctx.conflict_file}')
        conflict_states = ''.join(c[-1] for c in columns)
        if conflict_states != states:
            print(f'NOTE: {ctx.conflict_file} {window} states = {conflict_states}, '
                  f'microstate states = {states}; conflict L2 uses {conflict_states}.')
        conflict[f'conflict_{window}_l2'] = np.sqrt(
            np.square(conflict[columns].to_numpy(dtype=float)).sum(axis=1))

    # Keep a one-subject-one-row EEG feature table, without clinical columns.
    # Outer joins: Rest-based features are missing for sub-P001, but its
    # conflict and Simon duration features are kept.
    eeg = reorg[['subject_id'] + EEG_BASE].merge(
        conflict[['subject_id'] + EEG_CONFLICT], on='subject_id',
        how='outer', validate='one_to_one'
    ).merge(
        micro[['subject_id'] + EEG_DURATION], on='subject_id',
        how='left', validate='one_to_one')
    eeg.to_csv(ctx.csv('b6_eeg_features'), index=False)

    data = eeg.merge(
        subjects[['subject_id', 'group', 'age', 'gender', 'moca', 'updrs', 'trials_simon']],
        on='subject_id', how='inner', validate='one_to_one'
    ).merge(
        behavior[['subject_id', 'simon_err_inc']],
        on='subject_id', how='left', validate='one_to_one')

    data = apply_simon_qc(data, ctx)

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
        if len(pd_values) and len(hc_values):
            mw = mannwhitneyu(pd_values, hc_values, alternative='two-sided')
            mw_u, mw_p = float(mw.statistic), float(mw.pvalue)
        else:
            mw_u = mw_p = np.nan
        group_rows.append({'analysis':'PD_vs_Control', 'feature':feature,
                           'family':FAMILY, **output,
                           'mw_u': mw_u, 'mw_p': mw_p,
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
    gr.to_csv(ctx.csv('b6_group_stats'), index=False)
    cr.to_csv(ctx.csv('b6_correlations'), index=False)
    rq3.to_csv(ctx.csv('rq3_stats'), index=False)
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
    fig,ax=plt.subplots(figsize=(10.8, max(4.0, 0.47*len(EEG_FEATURES))))
    im=ax.imshow(heat,cmap='RdBu_r',vmin=-0.4,vmax=0.4,aspect='auto')
    ax.set_xticks(np.arange(len(keys)), ['MoCA\nAll','MoCA\nPD','Simon error\nAll','Simon error\nPD','UPDRS\nPD'])
    ax.set_yticks(np.arange(len(EEG_FEATURES)),EEG_FEATURES)
    for i in range(len(EEG_FEATURES)):
        for j in range(len(keys)):
            if np.isfinite(heat[i,j]): ax.text(j,i,f'{heat[i,j]:+.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,label='Age-adjusted partial Spearman rho')
    ax.set_title(f'B6: EEG microstate features and clinical/behavioral scores [{ctx.name}, K={len(states)}]')
    fig.tight_layout()
    fig.savefig(ctx.fig('bohyeon_b6_correlations'),dpi=220)
    plt.close(fig)
    save_run_config(ctx, states, 'B6_clinical_associations')
    print('\nOutputs saved in', RESULTS_DIR, 'and', FIGURES_DIR)


if __name__=='__main__':
    main()