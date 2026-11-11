
from fractions import Fraction as F
from itertools import product

def q(x):return x if isinstance(x,F) else F(str(x))
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(points):
 s=sorted(set(tuple(q(v) for v in p) for p in points))
 if len(s)<=1:return tuple(s)
 low=[];up=[]
 for p in s:
  while len(low)>=2 and cross(low[-2],low[-1],p)<=0:low.pop()
  low.append(p)
 for p in reversed(s):
  while len(up)>=2 and cross(up[-2],up[-1],p)<=0:up.pop()
  up.append(p)
 return tuple(low[:-1]+up[:-1])
def mat(A,x):return tuple(sum(q(A[i][j])*x[j] for j in range(2)) for i in range(2))
def plus(a,b):return tuple(a[k]+b[k] for k in range(2))
class Poly:
 def __init__(self,points):self.v=hull(points)
 @property
 def empty(self):return not bool(self.v)
 def affine(self,A):return Poly(mat(A,p) for p in self.v)
 def minkowski(self,other):return Poly(plus(p,z) for p in self.v for z in other.v)
 def bounds(self):return tuple((min(p[k] for p in self.v),max(p[k] for p in self.v)) for k in range(2))
 def contains_point(self,p):
  if self.empty:return False
  if len(self.v)==1:return tuple(p)==self.v[0]
  if len(self.v)==2:
   a,b=self.v
   return cross(a,b,p)==0 and all(min(a[k],b[k])<=p[k]<=max(a[k],b[k]) for k in range(2))
  return all(cross(self.v[i],self.v[(i+1)%len(self.v)],p)>=0 for i in range(len(self.v)))
 def contains_poly(self,other):return all(self.contains_point(p) for p in other.v)
 def clip_box(self,lo=(-1,-1),hi=(1,1)):
  ps=list(self.v)
  for n,c in [((q(1),q(0)),q(hi[0])),((q(-1),q(0)),-q(lo[0])),
              ((q(0),q(1)),q(hi[1])),((q(0),q(-1)),-q(lo[1]))]:
   if not ps:break
   dst=[]
   for a,b in zip(ps,ps[1:]+ps[:1]):
    fa=sum(n[k]*a[k] for k in range(2))-c
    fb=sum(n[k]*b[k] for k in range(2))-c
    ina=fa<=0;inb=fb<=0
    if ina:dst.append(a)
    if ina!=inb:
     frac=fa/(fa-fb)
     dst.append(tuple(a[k]+frac*(b[k]-a[k]) for k in range(2)))
   ps=dst
  return Poly(ps)
 def extrema(self,a):
  values=[sum(q(a[k])*p[k] for k in range(2)) for p in self.v]
  return min(values),max(values)
 def encode(self):return [[str(x) for x in p] for p in self.v]

def box2(x1lo,x1hi,x2lo,x2hi):
 return Poly(product((q(x1lo),q(x1hi)),(q(x2lo),q(x2hi))))

def reach(P,A,B,U,Wi):return P.affine(A).minkowski(U.affine(B)).minkowski(Wi)

def residual(P,Ai,Aj,Bi,Bj,U,Wi):
 DA=tuple(tuple(q(Ai[i][j])-q(Aj[i][j]) for j in range(2)) for i in range(2))
 DB=tuple(tuple(q(Bi[i][j])-q(Bj[i][j]) for j in range(2)) for i in range(2))
 return P.affine(DA).minkowski(U.affine(DB)).minkowski(Wi)
