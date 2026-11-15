"""Full offline pipeline smoke test with synthetic *reference* files only.
Not scientific validation; reports must not be interpreted as real Custo results."""
from pathlib import Path
import sys,struct,zipfile,io,shutil
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).parent))
import colab_pipeline as p
ROOT=Path(__file__).parent;WORK=ROOT/'_synthetic_e2e';shutil.rmtree(WORK,ignore_errors=True)
p.BASE=WORK;p.SOURCE=WORK/'source';p.RESULT=WORK/'results';p.FIG=WORK/'figures'
for q in [p.SOURCE,p.RESULT,p.FIG]:q.mkdir(parents=True,exist_ok=True)
d=np.load(ROOT/'tests/fixtures/reference_maps_rounded_1e-4.npz')
ch=pd.read_csv(ROOT/'tests/fixtures/k5_chan_coords_from_original_fif.csv')

def tag(kind,typ,payload):return struct.pack('>iiii',kind,typ,len(payload),0)+payload

def make_fif(filename,matrix,states):
    buf=bytearray()
    for _,r in ch.iterrows():
        payload=struct.pack('>iii',1,1,2)+struct.pack('>ff',1,1)+struct.pack('>i',1)+struct.pack('>12f',r.x_m,r.y_m,r.z_m,*([0]*9))+struct.pack('>ii',107,0)+r.channel.encode().ljust(16,b'\0')
        buf+=tag(203,30,payload)
    k,n=matrix.shape
    buf+=tag(3607,0x40000005,matrix.astype('>f8').tobytes()+struct.pack('>iii',n,k,2))
    buf+=tag(3502,10,':'.join(states).encode())
    (p.SOURCE/filename).write_bytes(buf)
for k in [4,5,6,7]:
    make_fif(f'templates_K{k}.fif',d[f'K{k}'],list('ABCDE'[:k]) if k==5 else (['A','B','C','Central'] if k==4 else list('ABCDEFG'[:k])))
make_fif('templates_K5_rest.fif',d['K5_rest'],list('ABCDE'))
make_fif('templates_K5_simon.fif',d['K5_simon'],list('ABCDE'))
if not (ROOT/'_parser_tests/synthetic_published_E.set').exists():
    import test_parsers  # generate the synthetic EEGLAB source fixture
shutil.copy(ROOT/'_parser_tests/synthetic_published_E.set',p.SOURCE/'Hwang 2025.set')
p.REFERENCE_FILENAMES=['Hwang 2025.set']
result=p.main()
assert len(result)==1 and (p.RESULT/'EXTERNAL_E_SHARED_REST_SIMON_COMPARISON.csv').exists()
print('PASS: offline synthetic end-to-end pipeline; full matrices+CSV+plots generated')
