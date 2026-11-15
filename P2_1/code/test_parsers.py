from pathlib import Path
import sys,struct,zipfile,io
import numpy as np
import pandas as pd
from scipy.io import savemat
sys.path.insert(0,str(Path(__file__).parent))
from colab_pipeline import load_team_modkmeans_fif,published_montage_and_solutions,analyse_one_reference
ROOT=Path(__file__).parent;P=ROOT/'_parser_tests';P.mkdir(exist_ok=True)
npz=np.load(ROOT/'tests/fixtures/reference_maps_rounded_1e-4.npz')
maps=npz['K5']; coords=pd.read_csv(ROOT/'tests/fixtures/k5_chan_coords_from_original_fif.csv')

def tag(kind,typ,p):return struct.pack('>iiii',kind,typ,len(p),0)+p
buf=bytearray()
for _,row in coords.iterrows():
    chname=row.channel.encode().ljust(16,b'\0')
    payload=(struct.pack('>iii',1,1,2)+struct.pack('>ff',1,1)+struct.pack('>i',1)
         +struct.pack('>12f',*[row.x_m,row.y_m,row.z_m]+[0]*9)
         +struct.pack('>ii',107,0)+chname)
    assert len(payload)==96
    buf+=tag(203,30,payload)
buf+=tag(3607,0x40000005,maps.astype('>f8').tobytes()+struct.pack('>iii',64,5,2))
buf+=tag(3502,10,b'A:B:C:D:E')
fif=P/'test_MNE_matrix.fif';fif.write_bytes(buf)
x,names,ch,xyz=load_team_modkmeans_fif(fif)
np.testing.assert_allclose(x,maps)
assert names==list('ABCDE') and ch==coords.channel.tolist()
print('PASS: synthetic FIF parser, 5x64 channel/order/positions/labels verified')

# Simulated EEGLAB nested msinfo with fourth map explicitly labeled E.
pub=pd.read_csv(ROOT/'tests/fixtures/Hwang2025_E_shared_electrodes_55.csv')
labels=pub.channel.tolist();n=len(labels)
refmaps=np.vstack([np.linspace(-1,1,n),np.cos(np.linspace(0,3,n)),np.sin(np.linspace(0,3,n)),pub.Hwang_2025_E.to_numpy()])
EEG={'chanlocs':np.array([{'labels':q} for q in labels],dtype=object),
     'msinfo':{'MSMaps':np.array([{'Maps':refmaps,'Labels':np.array(['A','B','C','E'],dtype=object)}],dtype=object)}}
path=P/'synthetic_published_E.set';savemat(path,{'EEG':EEG})
try:
    studies=published_montage_and_solutions(path)
    assert len(studies)==1
    assert studies[0]['labels']==['A','B','C','E']
    assert studies[0]['maps'].shape==(4,n)
    print('PASS: synthetic EEGLAB .set parser, explicit A-B-C-E labels preserved')
except Exception as exc:
    print('FAIL EEGLAB synthetic format',type(exc).__name__,str(exc))
    raise
