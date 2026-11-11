"""Small-dimensional EXACT rational polytopes, no floating-point optimisation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from typing import Iterable


def frac(x):
    if isinstance(x,F):return x
    return F(str(x))


def dot(a,b):
    if len(a)!=len(b):raise ValueError('Dimension mismatch in scalar product')
    return sum((x*y for x,y in zip(a,b)), F(0))


def gaussian_solve(A,b):
    n=len(b)
    aug=[[frac(A[i][j]) for j in range(n)]+[frac(b[i])] for i in range(n)]
    for c in range(n):
        k=next((r for r in range(c,n) if aug[r][c]!=0),None)
        if k is None:return None
        aug[c],aug[k]=aug[k],aug[c]
        pivot=aug[c][c]; aug[c]=[v/pivot for v in aug[c]]
        for r in range(n):
            if r!=c:
                z=aug[r][c]
                if z:aug[r]=[aug[r][j]-z*aug[c][j] for j in range(n+1)]
    return tuple(aug[i][n] for i in range(n))


@dataclass(frozen=True)
class Polytope:
    A:tuple
    b:tuple
    dim:int

    @classmethod
    def make(cls,A,b,box_lo,box_hi):
        lo=tuple(map(frac,box_lo)); hi=tuple(map(frac,box_hi))
        n=len(lo)
        if n==0 or len(hi)!=n or any(lo[i]>hi[i] for i in range(n)):
            raise ValueError('Nonempty dimension and valid finite enclosure required')
        aa=[tuple(map(frac,row)) for row in A]
        bb=[frac(x) for x in b]
        if len(aa)!=len(bb) or any(len(row)!=n for row in aa):raise ValueError('Invalid H-polytope')
        for j in range(n):
            unit=tuple(F(int(k==j)) for k in range(n))
            aa += [unit,tuple(-v for v in unit)]
            bb += [hi[j],-lo[j]]
        obj=cls(tuple(aa),tuple(bb),n)
        if not obj.vertices():raise ValueError('Empty polytope')
        return obj

    @classmethod
    def box(cls,lo,hi):
        return cls.make([],[],lo,hi)

    def contains(self,x):
        x=tuple(map(frac,x))
        if len(x)!=self.dim:return False
        return all(dot(row,x)<=bound for row,bound in zip(self.A,self.b))

    def vertices(self):
        out=set()
        for idx in combinations(range(len(self.b)),self.dim):
            x=gaussian_solve([self.A[j] for j in idx],[self.b[j] for j in idx])
            if x is not None and self.contains(x):out.add(x)
        return tuple(sorted(out))

    def intersects(self,other):
        if self.dim!=other.dim:raise ValueError('Dimension mismatch')
        A=self.A+other.A; b=self.b+other.b
        for idx in combinations(range(len(b)),self.dim):
            x=gaussian_solve([A[j] for j in idx],[b[j] for j in idx])
            if x is not None and all(dot(r,x)<=z for r,z in zip(A,b)):
                return True
        return False

    def shifted(self,delta):
        delta=tuple(map(frac,delta))
        if len(delta)!=self.dim:raise ValueError('Translation dimension')
        return Polytope(self.A,tuple(v+dot(row,delta) for row,v in zip(self.A,self.b)),self.dim)

    def public(self):
        return {'dim':self.dim,'H_rows':len(self.b),'vertices':[[str(z) for z in p] for p in self.vertices()]}


def components(polytopes):
    if not polytopes:raise ValueError('Nonempty family required')
    if any(not p.vertices() for p in polytopes.values()):raise ValueError('Nonempty compact polytopes required')
    dims={p.dim for p in polytopes.values()}
    if len(dims)!=1:raise ValueError('Mixed dimensions')
    names=list(polytopes)
    adj={k:set() for k in names};pairs=[]
    for a,b in combinations(names,2):
        intersect=polytopes[a].intersects(polytopes[b]);pairs.append((a,b,intersect))
        if intersect:adj[a].add(b);adj[b].add(a)
    unseen=set(names);parts=[]
    while unseen:
        start=sorted(unseen)[0];queue=[start];group=set()
        while queue:
            k=queue.pop()
            if k in group:continue
            group.add(k);queue.extend(adj[k]-group)
        unseen-=group;parts.append(tuple(sorted(group)))
    return {'connected_components':len(parts),
            'groups':sorted(parts),
            'intersecting_pairs':sum(bool(x[2]) for x in pairs),
            'disjoint_pairs':sum(not x[2] for x in pairs)}
