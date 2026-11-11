"""Independent floating-point finite first-exit facet LP verification.
"""
from scipy.optimize import linprog
import numpy as np

def solve(source,representative,X0,X,U,H):
    n=len(X.lo);m=len(U.lo);results=[]
    if H<1:raise ValueError('H>=1')
    for t in range(H):
        N=(t+1)*n+(t+1)*m+(t+1)*n
        def xi(k,p):return k*n+p
        def ui(k,p):return (t+1)*n+k*m+p
        def wi(k,p):return (t+1)*(n+m)+k*n+p
        eq=[];rhs=[];bounds=[]
        for k in range(t+1):
            for p in range(n):
                bounds.append((float(X0.lo[p]),float(X0.hi[p])) if k==0 else (float(X.lo[p]),float(X.hi[p])))
        for k in range(t+1):
            for p in range(m):bounds.append((float(U.lo[p]),float(U.hi[p])))
        for k in range(t+1):
            for p in range(n):bounds.append((float(source.W.lo[p]),float(source.W.hi[p])))
        for k in range(t):
            for p in range(n):
                row=np.zeros(N);row[xi(k+1,p)]=1
                for j in range(n):row[xi(k,j)]-=float(source.A[p][j])
                for j in range(m):row[ui(k,j)]-=float(source.B[p][j])
                row[wi(k,p)]-=1;eq.append(row);rhs.append(0)
        if eq:A_eq=np.array(eq);b_eq=np.array(rhs)
        else:A_eq=None;b_eq=None
        for r in range(n):
            for sign in (-1,+1):
                c=np.zeros(N)
                for j in range(n):c[xi(t,j)]-=sign*float(source.A[r][j]-representative.A[r][j])
                for j in range(m):c[ui(t,j)]-=sign*float(source.B[r][j]-representative.B[r][j])
                c[wi(t,r)]-=sign
                lp=linprog(c,A_eq=A_eq,b_eq=b_eq,bounds=bounds,method='highs')
                if lp.status==2:
                    results.append({'t':t,'facet':f'{sign:+d}{r}','status':'VACUOUS'})
                elif not lp.success:
                    raise RuntimeError(f'LP solver failed: {lp.message}')
                else:
                    h=representative.W.hi[r] if sign>0 else -representative.W.lo[r]
                    alpha=-lp.fun-float(h)
                    results.append({'t':t,'facet':f'{sign:+d}{r}','status':'NUMERICAL_PASS' if alpha<=1e-8 else 'NUMERICAL_VIOLATION','alpha':alpha})
    return results
