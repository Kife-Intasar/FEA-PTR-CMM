"""RQ1 reachable-prefix mode domination
"""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from random import Random
import json, time
from pathlib import Path

Q=lambda x:F(str(x))
A0=((Q('.5'),Q('.8')),(Q('.05'),Q('.4')))
X=(F(1),F(1)); C=(Q('.55'),Q('.16')); U=(Q('.1'),Q('.1'))
wj=(Q('.03'),Q('.005')); wi=(Q('.005'),Q('.005'))
Bplus=((F(1),Q('.1')),(Q('.1'),Q('.5')))
Bminus=((-F(1),Q('.1')),(Q('.1'),Q('.5')))

def addmat(A,D): return tuple(tuple(A[r][c]+D[r][c] for c in range(2)) for r in range(2))
def matvec(A,x): return tuple(sum(A[r][c]*x[c] for c in range(2)) for r in range(2))
def vecadd(*vs): return tuple(sum(v[i] for v in vs) for i in range(2))
def cert_invariance(m,C):
 A,B,w=m['A'],m['B'],m['w']
 return tuple(sum(abs(A[r][c])*C[c]+abs(B[r][c])*U[c] for c in range(2))+w[r] for r in range(2))
def containment(i,j,region):
 return tuple(i['w'][r]+sum(abs(i['A'][r][c]-j['A'][r][c])*region[c]+abs(i['B'][r][c]-j['B'][r][c])*U[c] for c in range(2))-j['w'][r] for r in range(2))

def build(k=6):
 modes=[]
 for sign,B in [(1,Bplus),(-1,Bminus)]:
  modes.append({'name':f'R{sign:+}', 'grp':sign,'A':A0,'B':B,'w':wj,'is_rep':True})
  for ind in range(k):
   lam=F(3,5)+F(ind, max(1,k-1))*F(2,5) if k>1 else F(1)
   D=((Q('.03')*lam,Q('.04')*lam),(F(0),F(0)))
   modes.append({'name':f'I{sign:+}_{ind+1}', 'grp':sign,'A':addmat(A0,D),'B':B,'w':wi,'is_rep':False})
 return modes

def solve2(B,y):
 det=B[0][0]*B[1][1]-B[0][1]*B[1][0]
 return ((y[0]*B[1][1]-y[1]*B[0][1])/det, (B[0][0]*y[1]-B[1][0]*y[0])/det)

def verify_constructive(modes):
 r=(Q('.04'),F(0)); u0=(Q('.1'),F(0))
 reps={m['grp']:m for m in modes if m['is_rep']}
 max_u=F(0); steps=0; biggest_violation=F(0)
 for m in modes:
  rep=reps[m['grp']]
  for w0s in product((-1,1),repeat=2):
   w0=tuple(F(w0s[a])*m['w'][a] for a in range(2))
   z=vecadd(matvec(m['B'],u0),w0)
   assert all(abs(z[a])<=C[a] for a in range(2))
   assert (z[0]>0)==(m['grp']==1) # distinguishes quotient classes
   v=solve2(rep['B'],tuple(r[a]-matvec(rep['A'],z)[a] for a in range(2)))
   assert all(abs(v[a])<=U[a] for a in range(2)),(m['name'],v)
   max_u=max(max_u,*(abs(v[a]) for a in range(2)))
   for w1s in product((-1,1),repeat=2):
    w1=tuple(F(w1s[a])*m['w'][a] for a in range(2))
    x2=vecadd(matvec(m['A'],z),matvec(m['B'],v),w1)
    assert all(abs(x2[a]-r[a])<=wj[a] for a in range(2)),(m['name'],x2)
    assert all(abs(x2[a])<=C[a] for a in range(2))
    steps+=1
 return {'vertex_cases':steps,'max_terminal_abs_u':float(max_u),'identified_quotient_modes':len(reps)}

def randomized_paths(modes, n=800):
 rng=Random(8741); reps={m['grp']:m for m in modes if m['is_rep']}
 checks=0
 for _ in range(n):
  i=rng.choice([m for m in modes if not m['is_rep']]);j=reps[i['grp']]
  x=(F(0),F(0))
  for t in range(5):
   u=tuple(F(rng.randint(-100,100),1000) for c in range(2))
   w=tuple(i['w'][c]*F(rng.randint(-100,100),100) for c in range(2))
   y=vecadd(matvec(i['A'],x),matvec(i['B'],u),w)
   inferred=tuple(y[c]-matvec(j['A'],x)[c]-matvec(j['B'],u)[c] for c in range(2))
   assert all(abs(inferred[c])<=j['w'][c] for c in range(2))
   assert all(abs(y[c])<=C[c] for c in range(2))
   x=y;checks+=1
 return checks

def lp_pair(i,j,h=4):
 try: import numpy as np;from scipy.optimize import linprog
 except ImportError: return {'available':False}
 n=2; st=2*(h+1); nu=2*(h+1); nw=2*(h+1)
 total=st+nu+nw
 def xi(t,r):return 2*t+r
 def ui(t,r):return st+2*t+r
 def wi(t,r):return st+nu+2*t+r
 eq=[];beq=[]
 for r in range(2):
  row=np.zeros(total);row[xi(0,r)]=1;eq.append(row);beq.append(0)
 for t in range(h):
  for r in range(2):
   row=np.zeros(total);row[xi(t+1,r)]=1
   for c in range(2): row[xi(t,c)]-=float(i['A'][r][c]);row[ui(t,c)]-=float(i['B'][r][c])
   row[wi(t,r)]-=1;eq.append(row);beq.append(0)
 bounds=[]
 for t in range(h+1):
  for r in range(2):bounds.append((-1,1))
 for t in range(h+1):
  for r in range(2):bounds.append((-float(U[r]),float(U[r])))
 for t in range(h+1):
  for r in range(2):bounds.append((-float(i['w'][r]),float(i['w'][r])))
 vals=[]
 for t in range(h+1):
  for r in range(2):
   for sgn in (-1,1):
    c=np.zeros(total)
    for a in range(2):c[xi(t,a)]-=sgn*float(i['A'][r][a]-j['A'][r][a]);c[ui(t,a)]-=sgn*float(i['B'][r][a]-j['B'][r][a])
    c[wi(t,r)]-=sgn
    res=linprog(c,A_eq=np.array(eq),b_eq=np.array(beq),bounds=bounds,method='highs')
    if not res.success:raise RuntimeError(res.message)
    vals.append(-res.fun-float(j['w'][r]))
 return {'available':True,'constraints_checked':len(vals),'max_margin':max(vals),'certified':max(vals)<1e-8}

def main():
 start=time.perf_counter()
 modes=build(k=6);reps={m['grp']:m for m in modes if m['is_rep']}
 inv={m['name']:cert_invariance(m,C) for m in modes}
 assert all(all(inv[m['name']][r]<=C[r] for r in range(2)) for m in modes)
 lr={m['name']:containment(m,reps[m['grp']],C) for m in modes if not m['is_rep']}
 gr={m['name']:containment(m,reps[m['grp']],X) for m in modes if not m['is_rep']}
 assert all(all(v<=0 for v in margins) for margins in lr.values())
 assert all(any(v>0 for v in margins) for margins in gr.values())
 corner=tuple(sum(reps[1]['A'][0][c] for c in range(2))-sum(abs(reps[1]['B'][0][c])*U[c] for c in range(2))-reps[1]['w'][0] for _ in range(1))
 assert corner[0]>1
 construct=verify_constructive(modes)
 rand=randomized_paths(modes)
 lps=[lp_pair(m,reps[m['grp']],h=3) for m in modes if not m['is_rep']]
 assert all(item.get('certified',False) for item in lps)
 opposite=[lp_pair(m,reps[-m['grp']],h=2) for m in modes if not m['is_rep']]
 assert all(not item['certified'] for item in opposite)
 q=6; scaling=[]
 for k in [1,2,3,4,5,6]:
  count=2*(k+1);q=2
  scaling.append({'modes':count,'global_supports':2**count-1,'reachable_supports':2**q-1,'enumeration_ratio':(2**count-1)/(2**q-1)})
 out={'system':{'X':['1','1'],'C':[str(t) for t in C],'U':[str(t) for t in U],'x0':['0','0'],'A_rep':[[str(v) for v in row] for row in A0]},'model_count':len(modes),'maximal_reachable_modes':2,'local_margin_max':float(max(max(z) for z in lr.values())),'global_margin_min':float(min(max(z) for z in gr.values())),'invariant_slack_min':float(min(C[r]-inv[m['name']][r] for m in modes for r in range(2))),'global_corner_min_next_first_coord':float(corner[0]),'constructive':construct,'randomized_exact_transition_checks':rand,'independent_lp':{'pairs':len(lps),'linear_programs':sum(z['constraints_checked'] for z in lps),'largest_margin':max(z['max_margin'] for z in lps),'all_passed':all(z['certified'] for z in lps),'negative_opposite_pairs':len(opposite),'negative_opposite_rejected':sum(not z['certified'] for z in opposite)},'support_scaling':scaling,'time_seconds':time.perf_counter()-start}
 Path(__file__).with_name('verification_results.json').write_text(json.dumps(out,indent=2))
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
