"""Independent Python verification of EEG microstate external map correspondence.

Run in Google Colab with Internet access. Requires only public GitHub FIF/.set files.
The output intentionally separates:  (1) exact common-named-electrode comparison,
(2) sensitivity check using standard_1005 spherical interpolation,
(3) optional subject-level bootstrap, if individual maps are supplied.

DO NOT cite (2) as a full native-coordinate Koenig/MSTemplateExplorer result.
DO NOT interpret spatial correlations as proof of identical neural generators.
"""
from __future__ import annotations

import io
import os
import re
import struct
import warnings
from pathlib import Path
from urllib.parse import quote

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.io import loadmat
from scipy.optimize import linear_sum_assignment
from mne.channels import make_standard_montage
from mne.channels.interpolation import _make_interpolation_matrix

try:
    import mat73
except ImportError:
    mat73 = None

BASE = Path('E_Microstate_External_Test')
BASE.mkdir(exist_ok=True)
SOURCE = BASE / 'source'
RESULT = BASE / 'results'
FIG = BASE / 'figures'
for p in [SOURCE, RESULT, FIG]:p.mkdir(exist_ok=True)

OUR_REPO = 'jimmytrn154/AI-based-Neural-Data-Science'
REF_REPO = 'ThomasKoenigBern/MS-Template-Explorer'
OUR_NAMES = ['A', 'B', 'C', 'D', 'E']

# Preselected map sources from public Koenig template database.
# Two K=5 sets (Qin), an E-labeled 4-state set (Hwang), Custo K=7,
# Zanesco datasets and 4-state examples are included to prevent label cherry-picking.
REFERENCE_FILENAMES = [
    'Custo 2017.set',
    'Zanesco 2020.set', 'Zanesco 2021a.set',
    'Hwang 2025.set',
    'Qin2022 Depression_Controls.set', 'Qin2022 Depression_Patients.set',
    'Zelenina_2022 Oxytocin.set',
    'Giordano 2018.set', 'Lassi2023.set', 'Milz_2016_Cognitive Style.set',
    'Schiller2021 Alcohol Intoxication.set', 'Toplutas 2023 DOC_HC.set',
]


def download(repo, path, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    url = f'https://raw.githubusercontent.com/{repo}/main/{quote(path, safe="/")}'
    response = requests.get(url, timeout=75)
    response.raise_for_status()
    dest.write_bytes(response.content)
    print('Downloaded', path, 'bytes=', len(response.content))
    return dest


def load_team_modkmeans_fif(path):
    """Read TEAM pycrostates ModKMeans .fif tags with stable MNE FIFF data layout.

    MNE FIFF: kind=203 channel information (96 bytes), kind=3607 PCA/cluster
    maps (dense float64 matrix), kind=3502 map labels. This is a *targeted*
    binary reader for this team's template files, not a generic FIF parser.
    """
    buf = Path(path).read_bytes()
    pos = 0
    channels = []
    xyz = []
    centers = None
    labels = None
    while pos + 16 <= len(buf):
        kind, type_code, size, next_pos = struct.unpack_from('>iiii', buf, pos)
        if size < 0 or pos + 16 + size > len(buf):
            raise ValueError(f'Malformed FIFF tag at {pos}')
        payload = memoryview(buf)[pos+16:pos+16+size]
        if kind == 203:   # ch_info struct
            channels.append(bytes(payload[80:96]).split(b'\0')[0].decode('utf-8'))
            xyz.append(struct.unpack_from('>fff', payload, 24))
        elif kind == 3607:  # dense double matrix of cluster centers
            ndim = struct.unpack_from('>i', payload, size-4)[0]
            if ndim != 2: raise ValueError(f'Unexpected center dims {ndim}')
            reverse_dims = struct.unpack_from('>ii', payload, size-12)
            dims = tuple(reversed(reverse_dims))
            rows, cols = dims
            expected_size = rows * cols * 8 + 12
            if expected_size != size: raise ValueError('Unexpected matrix size')
            centers = np.frombuffer(payload[:rows*cols*8], dtype='>f8').reshape(dims).astype(float)
        elif kind == 3502:
            labels = bytes(payload).decode('utf-8').split(':')
        if next_pos == -1: break
        pos = next_pos if next_pos > 0 else pos + 16 + size
    if centers is None or not channels or centers.shape[1] != len(channels):
        raise ValueError('Could not extract verified team template maps and channels')
    if labels is None or len(labels) != centers.shape[0]:
        raise ValueError('Cluster labels absent or wrong shape')
    xyz = np.asarray(xyz, dtype=float)
    if not np.isfinite(xyz).all() or np.any(np.linalg.norm(xyz,axis=1) < .001):
        raise ValueError('Channel digitization coordinates unusable')
    return centers, labels, channels, xyz


def upper_label(x):
    return str(x).strip().upper()


def normalize_label(obj):
    """Keep the publication's label, with a short final A-G label for summaries."""
    if isinstance(obj, (bytes, np.bytes_)):
        obj = bytes(obj).decode('utf-8', errors='replace')
    if isinstance(obj, (list, tuple, np.ndarray)):
        arr=np.asarray(obj, dtype=object).ravel().tolist()
        return normalize_label(arr[0]) if len(arr)==1 else str(arr)
    return str(obj).strip()


def published_letter(full):
    full = normalize_label(full)
    if re.fullmatch(r'[A-Z]', full):return full
    m = re.search(r'(?:CLASS|STATE|MICROSTATE|[_ -])([A-G])$',full,re.I)
    return m.group(1).upper() if m else full


def iterable_struct(obj):
    if obj is None:return []
    if isinstance(obj, np.ndarray):return obj.ravel().tolist()
    if isinstance(obj, (list,tuple)):return list(obj)
    return [obj]


def get_key(obj, key, default=None):
    if isinstance(obj, dict):return obj.get(key,default)
    return getattr(obj,key,default)


def published_montage_and_solutions(path):
    """Read EEGLAB .set file; support MATLAB v5 and v7.3 via scipy/mat73.

    Strict: use published MSMaps.Maps and MSMaps.Labels. Never infer label
    alphabetically from the map's ordinal position (Hwang fourth map is E).
    """
    path=Path(path)
    try:
        mat=loadmat(path, simplify_cells=True)
    except (NotImplementedError, ValueError, TypeError, OSError) as exc:
        if mat73 is None:
            raise RuntimeError('MATLAB v7.3 file needs pip install mat73') from exc
        mat=mat73.loadmat(str(path))
    eeg=mat.get('EEG')
    if eeg is None:raise ValueError('Missing EEG variable')
    if isinstance(eeg,(list,np.ndarray)) and len(eeg)==1:eeg=eeg[0]
    chanlocs=get_key(eeg,'chanlocs')
    channel_names=[]
    for ch in iterable_struct(chanlocs):
        channel_names.append(normalize_label(get_key(ch,'labels','')))
    if not channel_names or len(set(map(upper_label,channel_names)))!=len(channel_names):
        raise ValueError('Missing or duplicate valid EEG channel names')
    msinfo=get_key(eeg,'msinfo')
    mm=get_key(msinfo,'MSMaps')
    candidates=[]
    for sol in iterable_struct(mm):
        maps=get_key(sol,'Maps')
        labels=get_key(sol,'Labels')
        if maps is None or labels is None:continue
        try: maps=np.asarray(maps,dtype=float)
        except (ValueError,TypeError):continue
        if maps.ndim!=2 or maps.size==0:continue
        if maps.shape[1]==len(channel_names):pass
        elif maps.shape[0]==len(channel_names):maps=maps.T
        else:continue
        labels=[normalize_label(x) for x in iterable_struct(labels)]
        if maps.shape[0]!=len(labels) or not len(set(labels))==len(labels):
            # Some MATLAB v7.3 reader returns string chars in nested arrays;
            # do not silently assign A,B,C,D,E by position.
            continue
        if not (3<=maps.shape[0]<=12) or not np.isfinite(maps).all():continue
        candidates.append({'labels':labels,'maps':maps,'channels':channel_names,
                           'k':maps.shape[0], 'source_file':path.name})
    if not candidates:
        raise ValueError(f'No correctly labeled MSMaps in {path.name} (inspect MATLAB v7.3 structure)')
    return candidates


def abs_spatial(a,b):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if a.size!=b.size or a.size<3:return np.nan
    a=a-a.mean(); b=b-b.mean()
    d=np.linalg.norm(a)*np.linalg.norm(b)
    if not np.isfinite(d) or d==0:return np.nan
    return abs(float(a@b/d))


def matrix(a,b):
    return np.array([[abs_spatial(x,y) for y in b] for x in a])


def shared_label_matrix(ours, ours_ch, ref, ref_ch, min_electrodes=15):
    ours_lookup={upper_label(x):i for i,x in enumerate(ours_ch)}
    ref_lookup={upper_label(x):i for i,x in enumerate(ref_ch)}
    common=[c for c in ours_lookup if c in ref_lookup]
    if len(common)<min_electrodes:
        raise ValueError(f'Only {len(common)} common named channels; cannot make reliable named-only comparison')
    a=ours[:,[ours_lookup[c] for c in common]]
    b=ref[:,[ref_lookup[c] for c in common]]
    return matrix(a,b),common


STANDARD={str(k).upper():np.asarray(v) for k,v in
          make_standard_montage('standard_1005').get_positions()['ch_pos'].items()}


def standard_spline_sensitivity(ours, ours_ch, ref, ref_ch):
    """Approximate sensitivity only, NOT official coordinates-based result.

    Coordinates are standard_1005, not verified native reference digitizations.
    The comparison uses all recognizable channels of each template, not only
    the overlapping labels. It is designed to assess montage-induced changes.
    """
    ref_keep=[i for i,ch in enumerate(ref_ch) if upper_label(ch) in STANDARD]
    our_keep=[i for i,ch in enumerate(ours_ch) if upper_label(ch) in STANDARD]
    if len(ref_keep)<20 or len(our_keep)<20:
        raise ValueError('Insufficient standard_1005-identifiable coordinates')
    source=np.array([STANDARD[upper_label(ref_ch[i])] for i in ref_keep])
    target=np.array([STANDARD[upper_label(ours_ch[i])] for i in our_keep])
    W=_make_interpolation_matrix(source,target,alpha=1e-5)
    interp=ref[:,ref_keep]@W.T
    return matrix(ours[:,our_keep],interp),len(ref_keep),len(our_keep)


def one_to_one_assignment(sim, our_labels, external_labels):
    # Hungarian assignment over all candidate maps; the E-only best match
    # can differ from the globally optimal one-to-one solution.
    i,j=linear_sum_assignment(-sim)
    return [{'our_state':our_labels[a],'reference_state':external_labels[b],
             'abs_r':float(sim[a,b])} for a,b in zip(i,j)]


def analyse_one_reference(ours, our_labels, our_channels, source):
    refs=published_montage_and_solutions(source)
    all_rows=[]
    for item in refs:
        sim,shared=shared_label_matrix(ours,our_channels,item['maps'],item['channels'])
        full_labels=item['labels']
        short_labels=[published_letter(x) for x in full_labels]
        name=Path(source).stem
        outname=re.sub(r'[^a-zA-Z0-9_-]+','_',name)
        prefix=f'{outname}_K{item["k"]}'
        pd.DataFrame(sim,index=our_labels,columns=full_labels).to_csv(RESULT/f'{prefix}_absr.csv',float_format='%.8f')
        pd.DataFrame(100*sim**2,index=our_labels,columns=full_labels).to_csv(RESULT/f'{prefix}_shared_variance_pct.csv',float_format='%.4f')
        i=our_labels.index('E')
        order=np.argsort(-sim[i]); j0,j1=order[0],order[1]
        j_pub_e=[j for j,n in enumerate(short_labels) if n=='E']
        row={'file':Path(source).name,'reference_k':item['k'],
             'common_channels':len(shared),
             'best_reference_label':full_labels[j0],
             'best_reference_letter':short_labels[j0],
             'best_abs_r':sim[i,j0],
             'runner_up_reference_label':full_labels[j1],
             'runner_up_abs_r':sim[i,j1],
             'published_E_exists':bool(j_pub_e),
             'published_E_abs_r':sim[i,j_pub_e[0]] if j_pub_e else np.nan,
             'published_E_shared_variance_pct':100*sim[i,j_pub_e[0]]**2 if j_pub_e else np.nan,
             'E_margin_vs_best_competitor_sv_pct':(100*(sim[i,j_pub_e[0]]**2-
                max(sim[i,j]**2 for j in range(len(full_labels)) if j!=j_pub_e[0]))
                if j_pub_e else np.nan),
             'interpolation_method':'standard_1005_approx_diagnostic',
             'native_coordinate_matching_completed':False}
        try:
            sm,n_src,n_dst=standard_spline_sensitivity(ours,our_channels,item['maps'],item['channels'])
            row.update(approx_interp_best_letter=short_labels[np.argmax(sm[i])],
                       approx_interp_best_abs_r=float(np.max(sm[i])),
                       approx_interp_E_abs_r=(float(sm[i,j_pub_e[0]]) if j_pub_e else np.nan),
                       interp_ref_channel_count=n_src)
            pd.DataFrame(sm,index=our_labels,columns=full_labels).to_csv(RESULT/f'{prefix}_standard_spline_DIAGNOSTIC.csv',float_format='%.8f')
        except Exception as e:
            row['interp_error']=str(e)
        assignment=one_to_one_assignment(sim,our_labels,full_labels)
        row['hungarian_E_match']=next((x['reference_state'] for x in assignment if x['our_state']=='E'),None)
        all_rows.append(row)
        fig,ax=plt.subplots(figsize=(max(6,0.9*len(full_labels)+2),4.4))
        img=ax.imshow(100*sim**2,vmin=0,vmax=100)
        ax.set_xticks(range(len(full_labels)),short_labels)
        ax.set_yticks(range(len(our_labels)),our_labels)
        ax.set(xlabel=f'Published maps ({name})',ylabel='Our maps',title='Shared variance (%) • named-channel comparison')
        for rr in range(len(our_labels)):
            for cc in range(len(full_labels)):
                ax.text(cc,rr,f'{100*sim[rr,cc]**2:.0f}',ha='center',va='center',fontsize=8)
        fig.colorbar(img,ax=ax,label='Shared variance (%)');fig.tight_layout()
        fig.savefig(FIG/f'{prefix}_heatmap.png',dpi=160);plt.close(fig)
    return all_rows


def bootstrap_if_subject_maps_available(group_maps, names, channels, refs):
    """Optional exploratory bootstrap from individual *already-aligned* maps.

    Not equivalent to recomputing participant-level GFP clustering + second
    level k-means on every resample. Report as approximate bootstrap only.
    File: subject_maps.npz with maps=(subjects,5,64), groups=optional.
    """
    path=SOURCE/'subject_maps.npz'
    if not path.exists():
        print('SKIP participant bootstrap: subject_maps.npz absent in public repository/uploads.')
        return None
    z=np.load(path,allow_pickle=False)
    x=np.asarray(z['maps'],float)
    if x.ndim!=3 or x.shape[1:]!=group_maps.shape:
        raise ValueError('subject_maps.npz maps must be (n_subjects, 5, 64) aligned to A–E')
    if not np.isfinite(x).all():raise ValueError('nonfinite subject maps')
    # Align each subject's entire map set using Hungarian; sign flip when needed.
    aligned=[]
    for subject in x:
        sim=matrix(group_maps,subject)
        rr,cc=linear_sum_assignment(-sim)
        target=np.empty_like(subject)
        for i,j in zip(rr,cc):
            v=subject[j].copy()
            if np.corrcoef(v,group_maps[i])[0,1]<0:v=-v
            target[i]=v/np.linalg.norm(v)
        aligned.append(target)
    aligned=np.asarray(aligned)
    rng=np.random.default_rng(20261110)
    n=len(aligned); reps=1000
    groups = np.asarray(z['groups']).astype(str) if 'groups' in z.files else None
    strata = [np.flatnonzero(groups==g) for g in sorted(set(groups))] if groups is not None else [np.arange(n)]
    summary=[]
    for ref in refs:
        refmap=ref['maps'];refch=ref['channels'];labels=ref['labels']
        jE=[j for j,s in enumerate(labels) if published_letter(s)=='E']
        if not jE:continue
        bootstrap_matches=[];sv=[]
        for _ in range(reps):
            picks=np.concatenate([rng.choice(idx,size=len(idx),replace=True) for idx in strata])
            prototype=np.mean(aligned[picks,4,:],axis=0)
            sim,_=shared_label_matrix(prototype[None],channels,refmap,refch,min_electrodes=15)
            bootstrap_matches.append(int(np.argmax(sim[0])==jE[0]))
            sv.append(float(sim[0,jE[0]]**2))
        summary.append({'ref':ref['source_file'],'approx_bootstrap_reps':reps,'n_subjects':n,
                        'E_is_best_fraction':np.mean(bootstrap_matches),
                        'E_shared_variance_CI_low_pct':100*np.quantile(sv,.025),
                        'E_shared_variance_CI_high_pct':100*np.quantile(sv,.975)})
    if summary:pd.DataFrame(summary).to_csv(RESULT/'OPTIONAL_approx_subject_bootstrap.csv',index=False)
    return summary


def main():
    paths={}
    for k in [4,5,6,7]:
        paths[k]=download(OUR_REPO,f'results/templates/templates_K{k}.fif',SOURCE/f'templates_K{k}.fif')
    rest=download(OUR_REPO,'results/templates/templates_K5_rest.fif',SOURCE/'templates_K5_rest.fif')
    simon=download(OUR_REPO,'results/templates/templates_K5_simon.fif',SOURCE/'templates_K5_simon.fif')
    own={}
    for k,path in paths.items():
        own[k]=load_team_modkmeans_fif(path)
        print(f'OUR K{k}:',own[k][0].shape,'channels',len(own[k][2]))
    maps5,names,channels,xyz=own[5]
    if names!=OUR_NAMES:raise ValueError('Expected canonical team labels A–E; inspect original file')
    pd.DataFrame(xyz,columns=['x_m','y_m','z_m']).assign(channel=channels).to_csv(
        RESULT/'our_K5_original_channel_positions.csv',index=False)
    internal=matrix(maps5,maps5)
    np.testing.assert_allclose(internal[4,2],0.336653,atol=2e-5)
    pd.DataFrame(internal,index=names,columns=names).to_csv(RESULT/'our_internal_K5_absr.csv')
    compareK4=matrix(own[4][0],maps5)
    pd.DataFrame(compareK4,index=own[4][1],columns=names).to_csv(RESULT/'our_K4_vs_K5_absr.csv')
    for k in [6,7]:
        pd.DataFrame(matrix(maps5,own[k][0]),index=names,columns=own[k][1]).to_csv(
            RESULT/f'our_K5_vs_K{k}_absr.csv')
    rmap=load_team_modkmeans_fif(rest)[0]
    smap=load_team_modkmeans_fif(simon)[0]
    print('Internal Rest E vs Simon E:',matrix(rmap,smap)[4,4])
    failures=[]; rows=[];bootstrap_refs=[];condition_rows=[]
    for filename in REFERENCE_FILENAMES:
        p=SOURCE/filename
        try:
            download(REF_REPO,'Template Maps/'+filename,p)
            studies=published_montage_and_solutions(p)
            bootstrap_refs.extend(studies)
            rows.extend(analyse_one_reference(maps5,names,channels,p))
            for item in studies:
                for condition, maps_here in [('shared',maps5),('rest_only',rmap),('simon_only',smap)]:
                    sim,shared=shared_label_matrix(maps_here,channels,item['maps'],item['channels'])
                    ii=names.index('E'); candidate_names=[published_letter(q) for q in item['labels']]
                    jE=[j for j,s in enumerate(candidate_names) if s=='E']
                    condition_rows.append({
                        'reference_file':filename,'reference_k':item['k'],'our_map_set':condition,
                        'n_common_channels':len(shared),
                        'best_reference_label':item['labels'][int(np.argmax(sim[ii]))],
                        'best_abs_r':float(np.max(sim[ii])),
                        'published_E_abs_r':float(sim[ii,jE[0]]) if jE else np.nan})
        except Exception as exc:
            err=str(exc)
            print('FAILED external reference:',filename,err)
            failures.append({'filename':filename,'error':err})
    if rows:
        out=pd.DataFrame(rows)
        out.to_csv(RESULT/'EXTERNAL_MATCHING_ALL_STUDIES.csv',index=False,float_format='%.7f')
        print('\nOur E versus all published maps (exact named-channel baseline)')
        print(out[['file','reference_k','common_channels','best_reference_letter','best_abs_r',
                   'published_E_exists','published_E_abs_r','approx_interp_best_letter']].to_string(index=False))
    if condition_rows:
        pd.DataFrame(condition_rows).to_csv(RESULT/'EXTERNAL_E_SHARED_REST_SIMON_COMPARISON.csv',index=False)
    if failures:
        pd.DataFrame(failures).to_csv(RESULT/'REFERENCES_NOT_PARSED.csv',index=False)
        print('WARNING: Failed to parse some sources; do not claim the reference survey is exhaustive.')
    if not rows:raise RuntimeError('No published map was successfully parsed: inspect downloaded .set files.')
    bootstrap_if_subject_maps_available(maps5,names,channels,bootstrap_refs)
    print('\nSUCCESS. Output directory:',RESULT.resolve())
    print('Important: standard_1005 interpolation is a sensitivity check, not verified native-coordinate matching.')
    return out if rows else None

if __name__=='__main__':main()
