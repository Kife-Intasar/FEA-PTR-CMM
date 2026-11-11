"""Exact rational interval operations.
"""
from dataclasses import dataclass
from fractions import Fraction as Q

def f(v):
    return v if isinstance(v,Q) else Q(str(v))

@dataclass(frozen=True)
class Box:
    lo: tuple[Q,...]
    hi: tuple[Q,...]
    def __post_init__(self):
        if not (len(self.lo)==len(self.hi) and all(a<=b for a,b in zip(self.lo,self.hi))):
            raise ValueError('Invalid box')
    @classmethod
    def make(cls,lo,hi):return cls(tuple(map(f,lo)),tuple(map(f,hi)))
    @classmethod
    def symmetric(cls,width):
        w=tuple(map(f,width));return cls(tuple(-v for v in w),w)
    @classmethod
    def singleton(cls,values):
        a=tuple(map(f,values));return cls(a,a)
    def intersect(self,other):
        a=[max(x,y) for x,y in zip(self.lo,other.lo)]
        b=[min(x,y) for x,y in zip(self.hi,other.hi)]
        if any(x>y for x,y in zip(a,b)):return None
        return Box(tuple(a),tuple(b))
    def within(self,other):
        return len(self.lo)==len(other.lo) and all(a>=c and b<=d for a,b,c,d in zip(self.lo,self.hi,other.lo,other.hi))
    def export(self):return {'lo':list(map(str,self.lo)), 'hi':list(map(str,self.hi))}

def matrix(A):return tuple(tuple(map(f,r)) for r in A)

def linear(A,X):
    lo=[];hi=[]
    for row in A:
        lo.append(sum((z*(a if z>=0 else b) for z,a,b in zip(row,X.lo,X.hi)),Q(0)))
        hi.append(sum((z*(b if z>=0 else a) for z,a,b in zip(row,X.lo,X.hi)),Q(0)))
    return Box(tuple(lo),tuple(hi))

def add(*args):
    return Box(tuple(sum((box.lo[i] for box in args),Q(0)) for i in range(len(args[0].lo))),
               tuple(sum((box.hi[i] for box in args),Q(0)) for i in range(len(args[0].hi))))

def subtract(A,B):return tuple(tuple(a-b for a,b in zip(x,y)) for x,y in zip(A,B))

@dataclass(frozen=True)
class Mode:
    name: str
    A: tuple
    B: tuple
    W: Box
    @classmethod
    def make(cls,name,A,B,W):
        a=matrix(A);b=matrix(B)
        if len(a)!=len(W.lo) or any(len(row)!=len(a) for row in a) or len(b)!=len(a):raise ValueError('Dimensions')
        return cls(name,a,b,W)

def successor(mode,X,U):return add(linear(mode.A,X),linear(mode.B,U),mode.W)
def residual(source,rep,X,U):return add(linear(subtract(source.A,rep.A),X),linear(subtract(source.B,rep.B),U),source.W)

def finite_outer(source,representative,X0,X,U,H):
    reached=X0.intersect(X)
    if reached is None: return {'status':'PASS_VACUOUS','H':H,'steps':[],'safe_terminal':None}
    steps=[]
    for t in range(H):
        r=residual(source,representative,reached,U)
        good=r.within(representative.W)
        steps.append({'t':t,'safe_outer':reached.export(), 'residual_outer':r.export(), 'rational_pass':good})
        if not good:
            return {'status':'UNKNOWN_OUTER','H':H,'steps':steps,'safe_terminal':None}
        reached=successor(source,reached,U).intersect(X)
        if reached is None:
            return {'status':'PASS_EMPTY_TERMINAL','H':H,'steps':steps,'safe_terminal':None}
    return {'status':'PASS','H':H,'steps':steps,'safe_terminal':reached}

def continuation(source,rep,terminal,X,U,C):
    if terminal is None:
        return {'status':'CERTIFIED_EMPTY_TAIL','capture':True,'closure':True,'residual':True}
    if C is None:return {'status':'UNKNOWN_NO_GUARD','capture':False,'closure':False,'residual':False}
    cap=terminal.within(C) and C.within(X)
    succ=successor(source,C,U)
    intersected=succ.intersect(X)
    closed=intersected is None or intersected.within(C)
    resid=residual(source,rep,C,U).within(rep.W)
    return {'status':'CERTIFIED' if (cap and closed and resid) else 'UNKNOWN', 'capture':cap,'closure':closed,'residual':resid,
            'guard':C.export(),'successor_outer':succ.export()}
