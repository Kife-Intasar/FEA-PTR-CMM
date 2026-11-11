"""Task-specific FINITE-STATE instantiations.
Krook et al 2024: Definition 1 enforceable superblock strong/weak-until characterisation.
Ganty et al 2025: exact reachable-source simulation preorder on a finite LTS.
Ganty et al 2026: finite powerset-closure reachability-region semi-complete refinement.
"""
from dataclasses import dataclass
from itertools import product
from collections import deque

@dataclass(frozen=True)
class LTS:
    states: tuple
    actions: tuple
    edges: dict
    labels: tuple
    initial: tuple

    def __post_init__(self):
        assert self.states == tuple(range(len(self.states)))
        assert len(self.labels)==len(self.states)
        for s in self.states:
            for a in self.actions:
                v=self.edges.get((s,a),())
                assert v and all(x in self.states for x in v), (s,a,v)

    @property
    def allbits(self):return (1<<len(self.states))-1
    def mask(self,seq):return sum(1<<x for x in set(seq))
    def reachable(self):
        seen=set(self.initial);q=deque(self.initial)
        while q:
            s=q.popleft()
            for a in self.actions:
                for t in self.edges[s,a]:
                    if t not in seen:seen.add(t);q.append(t)
        return frozenset(seen)

def partition_by_labels(ts):
    groups={}
    for s in ts.states:groups.setdefault(ts.labels[s], []).append(s)
    return tuple(tuple(v) for k,v in sorted(groups.items(),key=lambda z:str(z[0])))

def controllable_until(ts,P,T, strong):
    result=set(T) if strong else set(P)|set(T)
    while True:
        new=set(T)|{s for s in P if any(set(ts.edges[s,a])<=result for a in ts.actions)}
        if new==result:return frozenset(result)
        result=new

def krook_rsbs_partition(ts):
    part=partition_by_labels(ts);split_count=0;iterations=0;formulae=0
    while True:
        iterations+=1;changed=False
        for bi,P in enumerate(part):
            other=[part[j] for j in range(len(part)) if j!=bi]
            signatures={s:[] for s in P}
            for k in range(1<<len(other)):
                T=set().union(*(other[j] for j in range(len(other)) if k>>j&1)) if k else set()
                for strong in (True,False):
                    winning=controllable_until(ts,set(P),T,strong)
                    formulae+=1
                    for s in P:signatures[s].append(s in winning)
            groups={}
            for s in P:groups.setdefault(tuple(signatures[s]),[]).append(s)
            if len(groups)>1:
                newblocks=tuple(tuple(g) for g in groups.values())
                part=part[:bi]+newblocks+part[bi+1:]
                split_count+=len(newblocks)-1;changed=True;break
        if not changed:break
        assert iterations<=len(ts.states)+1
    # Explicit verify Definition 1 for each pair and all P,U,T strong/weak superblocks
    for P in part:
        others=[part[k] for k,Q in enumerate(part) if Q!=P]
        for mask in range(1<<len(others)):
            T=set().union(*(others[k] for k in range(len(others)) if mask>>k&1)) if mask else set()
            for strong in (True,False):
                W=controllable_until(ts,set(P),T,strong)
                assert len({s in W for s in P})==1
    return {'partition':part,'blocks':len(part),'splits':split_count,'iterations':iterations,'formula_evaluations':formulae,
            'definition1_checked':True}

def reachable_simulation_2025(ts):
    sources=ts.reachable()
    relation={(p,q) for p in sources for q in ts.states if ts.labels[p]==ts.labels[q]}
    start=len(relation);loops=0
    while True:
        loops+=1
        refined=set()
        for p,q in relation:
            if all(all(any((p2,q2) in relation for q2 in ts.edges[q,a])
                           for p2 in ts.edges[p,a]) for a in ts.actions):
                refined.add((p,q))
        if refined==relation:break
        relation=refined
    assert all((s,s) in relation for s in sources)
    return {'relation':relation,'reachable_states':sources,'simulating_pairs':len(relation),
            'candidate_pairs':start,'iterations':loops,'identity_pairs':len(sources)}

def universal_pre(ts,mask):
    out=0
    for s in ts.states:
        if all(all(mask>>x&1 for x in ts.edges[s,a]) for a in ts.actions):out|=1<<s
    return out

def closure(part,mask):
    return sum(sum(1<<x for x in P) for P in part if any(mask>>x&1 for x in P))

def semi_complete_2026(ts,relevant=None):
    n=len(ts.states)
    assert n<=11,'exhaustive finite-state study restricted to n<=11'
    relevant=set(ts.states) if relevant is None else set(relevant)
    part=partition_by_labels(ts);initial=len(part);splits=0;rounds=0;checks=0
    rel_true=ts.mask(relevant)
    while True:
        rounds+=1
        rel_rounded=closure(part,rel_true)  
        violation=None
        for mask in range(1<<len(part)):
            abstract=sum(ts.mask(part[k]) for k in range(len(part)) if mask>>k&1)
            right=universal_pre(ts,abstract)&rel_rounded  
            left=closure(part,right)&rel_rounded          
            checks+=1
            if left!=right:
                violation=(mask,right,left,rel_rounded)
                break
        if violation is None:break
        _,right,left,rel_rounded=violation
        selected=None
        for j,P in enumerate(part):
            if not any(rel_rounded>>s&1 for s in P):continue
            a=tuple(s for s in P if right>>s&1)
            b=tuple(s for s in P if not(right>>s&1))
            if a and b:selected=j,a,b;break
        assert selected is not None,(part,violation)
        j,a,b=selected
        part=part[:j]+(a,b)+part[j+1:]
        splits+=1
        assert splits<=n-initial
    # Verify EXACT Definition 3.3 on every abstract value at final partition.
    rel_rounded=closure(part,rel_true)
    for mask in range(1<<len(part)):
        abstract=sum(ts.mask(part[k]) for k in range(len(part)) if mask>>k&1)
        right=universal_pre(ts,abstract)&rel_rounded
        left=closure(part,right)&rel_rounded
        assert left==right
    return {'partition':part,'blocks':len(part),'splits':splits,'iterations':rounds,
            'abstract_value_checks':checks,'initial_blocks':initial,
            'reachable_states':len(relevant),'rounded_relevant_states':rel_rounded.bit_count(),
            'paper_def_3_3_ideal_semi_completeness_verified':True}
def grid_mimo_fixture(which):
    from itertools import product
    points=list(product((-1,0,1),repeat=2)); n=len(points);exit_state=n
    index={p:i for i,p in enumerate(points)}
    if which=='coupled_robust':
        A=((0,1),(0,0));actions=((0,0),(1,0),(-1,0));noise=((0,0),(0,1))
    elif which=='coupled_sparse':
        A=((0,1),(0,0));actions=((0,0),);noise=((0,0),)
    elif which=='coupled_recurrent':
        A=((0,1),(1,0));actions=((0,0),);noise=((0,0),(0,1))
    else:raise ValueError(which)
    edges={}
    for s,p in enumerate(points):
        for u in actions:
            out=set()
            for w in noise:
                z=(A[0][0]*p[0]+A[0][1]*p[1]+u[0]+w[0],
                   A[1][0]*p[0]+A[1][1]*p[1]+u[1]+w[1])
                out.add(index.get(z,exit_state))
            edges[s,u]=tuple(sorted(out))
    for u in actions:edges[exit_state,u]=(exit_state,)
    labels=tuple('GOAL' if p==(0,0) else 'SAFE' for p in points)+('EXIT',)
    return LTS(tuple(range(n+1)),actions,edges,labels,(index[(1,0)],))

def unreachable_component_control():
    states=tuple(range(8));action=('tick',)
    successor={0:1,1:1,2:3,3:4,4:4,5:6,6:7,7:7}
    edges={(s,'tick'):(successor[s],) for s in states}
    labels=tuple('GOAL' if s==1 else 'EXIT' if s==7 else 'SAFE' for s in states)
    return LTS(states,action,edges,labels,(0,))

def v_interval_gridded_fixture(mode):
    from fractions import Fraction as F
    intervals=((F(-1),-F(1,5)),(-F(1,5),F(1,5)),(F(1,5),F(1)))
    cells=list(product(intervals,repeat=2));index={z:k for k,z in enumerate(cells)}
    exit_state=len(cells)
    actions=((F(-1),F(-1)),(F(-1),F(1)),(F(0),F(0)),(F(1),F(-1)),(F(1),F(1)))
    eps=F(1,10000)
    if mode=='source':
        A=((F(1),F(-1)),(F(1),F(-1)))
        W=((-eps,eps),(-eps,eps))
    elif mode=='representative':
        A=((F(41,40),-F(21,20)),(F(1),F(-1)))
        W=((-F(1,40)-eps,F(11,400)+F(41,40)*eps),(-eps,eps))
    else:raise ValueError(mode)
    def successor_box(cell,u):
        result=[]
        for k,row in enumerate(A):
            lower=sum(v*(cell[j][0] if v>=0 else cell[j][1]) for j,v in enumerate(row))
            upper=sum(v*(cell[j][1] if v>=0 else cell[j][0]) for j,v in enumerate(row))
            result.append((lower+u[k]/10+W[k][0],upper+u[k]/10+W[k][1]))
        return result
    def intersects(a,b):return a[0]<=b[1] and b[0]<=a[1]
    edges={}
    for k,cell in enumerate(cells):
        for u in actions:
            res=successor_box(cell,u)
            nxt={q for q,dest in enumerate(cells) if intersects(res[0],dest[0]) and intersects(res[1],dest[1])}
            if any(c[0]<-1 or c[1]>1 for c in res):nxt.add(exit_state)
            assert nxt
            edges[k,u]=tuple(sorted(nxt))
    for u in actions:edges[exit_state,u]=(exit_state,)
    labels=tuple('GOAL_CELL' if i==4 else 'SAFE_CELL' for i in range(9))+('EXIT',)
    initial=tuple(i for i,c in enumerate(cells) if c[0][0]<=1<=c[0][1] and c[1][0]<=0<=c[1][1])
    assert initial
    return LTS(tuple(range(10)),actions,edges,labels,initial)
