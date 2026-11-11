"""Production of uploaded fourteen-mode source."""
from .exactbox import Mode,Box,f

def fourteen():
    from pathlib import Path
    import sys
    base=Path(__file__).resolve().parents[1]/'source'
    sys.path.insert(0,str(base))
    from canonical_verify_rq1_reachable import build,C,U,X
    M=build(6)
    mods={d['name']:Mode.make(d['name'],d['A'],d['B'],Box.symmetric(d['w'])) for d in M}
    reps={d['grp']:d['name'] for d in M if d['is_rep']}
    pair=[(d['name'],reps[d['grp']]) for d in M if not d['is_rep']]
    return mods,pair,Box.singleton([0,0]),Box.symmetric(X),Box.symmetric(U),Box.symmetric(C)

def bridge():
    A_i=[['1/2','1/4'],['0','1/2']];A_j=[['11/20','1/4'],['0','1/2']]
    B=[['1','0'],['0','1']]
    return (Mode.make('B_i',A_i,B,Box.symmetric(['1/100','1/100'])),
            Mode.make('B_j',A_j,B,Box.symmetric(['1/10','1/10'])),
            Box.singleton(['1','0']),Box.symmetric(['2','2']),Box.symmetric(['1/100','1/100']),Box.symmetric(['1/4','1/4']))

def strict():
    A_i=[['-2','0'],['1/5','0']];A_j=[['-1','0'],['1/5','0']]
    B=[['1/100','1/100'],['1/100','-1/100']]
    return (Mode.make('S_i',A_i,B,Box.singleton([0,0])),
            Mode.make('S_j',A_j,B,Box.make(['-9/10','0'],['1/2','0'])),
            Box.singleton(['1/5','0']),Box.symmetric(['1','1']),Box.symmetric(['1','1']))
