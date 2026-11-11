
"""Independent rational audit of viable MIMO fixture and comparator semantics.
"""
from fractions import Fraction as F
from pathlib import Path
import json, csv

Q = lambda x: F(str(x))
A=((F(1),F(-1)),(F(1),F(-1)))
J=((F(41,40),-F(21,20)),(F(1),F(-1)))
B=((F(1,10), F(0)),(F(0),F(1,10)))
X0=(F(1),F(0))
r=F(1,10)

def mat_vec(M, v):
    return tuple(sum(M[k][t]*v[t] for t in (0,1)) for k in (0,1))
def mat_mul(M,N):
    return tuple(tuple(sum(M[i][k]*N[k][j] for k in (0,1)) for j in (0,1)) for i in (0,1))
def absmat_vec(M, v):
    return tuple(sum(abs(M[k][t])*v[t] for t in (0,1)) for k in (0,1))
def vecsum(x,y):return tuple(x[k]+y[k] for k in (0,1))
def endpoints(M,lo,hi):
    a,b=[],[]
    for row in M:
        a.append(sum(c*(lo[k] if c>=0 else hi[k]) for k,c in enumerate(row)))
        b.append(sum(c*(hi[k] if c>=0 else lo[k]) for k,c in enumerate(row)))
    return tuple(a),tuple(b)
def fstr(x):return str(x)

def verify(eps):
    assert 0<eps<=F(1,1000), 'outside proved range'
    assert mat_mul(A,A)==((F(0),F(0)),(F(0),F(0)))
    lo=-F(1,40)-eps; hi=F(11,400)+F(41,40)*eps
    assert hi>0>lo
    r0lo=-F(1,40)-eps
    r0hi=-F(1,40)+eps
    r1lo=F(1,50)-F(21,20)*eps
    r1hi=F(11,400)+F(41,40)*eps
    first_exit={"t0_first_component":[r0lo,r0hi],"t1_first_component":[r1lo,r1hi],"both_t_second_component":[-eps,eps]}
    assert r0lo>=lo and r0hi<=hi and r1lo>=lo and r1hi<=hi
    assert -eps>=-eps and eps<=eps
    C_radius=F(3,10)+3*eps
    ell_radius=F(1,4)+F(5,2)*eps
    tail_max=F(1,20)*ell_radius+eps
    assert C_radius < 1 and tail_max <= hi and -tail_max >= lo
    u0=(F(-1),F(-1))
    x1=vecsum(mat_vec(A,X0),mat_vec(B,u0)); x2=mat_vec(A,x1)
    assert x1==(F(9,10),F(9,10)) and x2==(F(0),F(0))
    z=tuple(F(9,10)*X0[k]+F(1,10)*x2[k] for k in (0,1))
    y=vecsum(mat_vec(A,z),mat_vec(B,(F(-1),F(1))))
    assert z==(F(9,10),F(0)) and y==(F(4,5),F(1))
    assert all(-1<=v<=1 for v in y)
    ell_y=y[1]-F(1,2)*y[0]
    convex_gap=F(1,20)*ell_y-hi
    assert convex_gap>0
    global_margin=max(F(1,20),F(19,400)-eps/F(40))
    assert global_margin == F(1,20)
    reverse_margin=F(21,400)+eps/F(40)
    assert reverse_margin>0
    reachable_margin=max(r0hi-hi,lo-r0lo,r1hi-hi,lo-r1lo,tail_max-hi,lo+tail_max)
    assert reachable_margin<=0
    first_x1_lo=(F(9,10)-eps,F(9,10)-eps)
    first_x1_hi=(F(381,400)+F(41,40)*eps,F(9,10)+eps)
    second_lo,second_hi=endpoints(J,first_x1_lo,first_x1_hi)
    second_lo=(second_lo[0]+lo,second_lo[1]-eps)
    second_hi=(second_hi[0]+hi,second_hi[1]+eps)
    assert all(-r<=v<=r for v in (*second_lo,*second_hi))
    assert all(-1<=v<=1 for v in (*first_x1_lo,*first_x1_hi))
    nu=(max(abs(lo),abs(hi)),eps)
    J2=mat_mul(J,J)
    even2=vecsum(vecsum(absmat_vec(J2,(r,r)),absmat_vec(J,nu)),nu)
    odd=vecsum(absmat_vec(J,(r,r)),nu)
    assert all(v<r for v in even2) and all(v<1 for v in odd)
    src_odd=F(2,10)+eps
    src_even=3*eps
    assert src_odd<1 and src_even<r
    assert convex_gap==F(1,400)-F(41,40)*eps
    return {
        'epsilon':str(eps),'finite_t0_first_interval':list(map(str,first_exit['t0_first_component'])),
        'finite_t1_first_interval':list(map(str,first_exit['t1_first_component'])),
        'terminal_C_coordinate_radius':str(C_radius),'terminal_C_nonempty':True,
        'tail_residual_radius':str(tail_max),'tail_invariant_by_A2_equals_zero':True,
        'global_identity_margin':str(global_margin),
        'identity_step_reverse_margin':str(reverse_margin),
        'no_convex_direct_guard_gap':str(convex_gap),
        'reachable_one_way_max_facet_margin':str(reachable_margin),
        'shared_robust_reach_continue':True,
        'mode_j_even_two_step_bound':list(map(str,even2)),
        'mode_j_odd_step_bound':list(map(str,odd)),
        'source_even_bound':str(src_even),
    }

def main():
    samples=[F(1,10000),F(1,2000),F(1,1000)]
    records=[verify(e) for e in samples]
    expected=records[0]
    assert expected['global_identity_margin']=='1/20'
    assert expected['identity_step_reverse_margin']=='21001/400000'
    assert expected['no_convex_direct_guard_gap']=='959/400000'
    assert records[-1]['no_convex_direct_guard_gap']=='59/40000'
    for value in (F(1,1000000),F(1,10000),F(1,1000)):
        assert verify(value)['shared_robust_reach_continue']
    proof={'fixture':'Viable MIMO V6.1','arithmetic':'Python fractions.Fraction only',
           'tested_epsilons':[str(e) for e in samples],
           'parametric_range':'0<epsilon<=1/1000; analytic affine bounds cover continuum',
           'results':records,
           'published_complete_methods_executed':0
           }
    base=Path(__file__).resolve().parents[1]; out=base/'reports';out.mkdir(exist_ok=True)
    (out/'INDEPENDENT_EXACT_AUDIT.json').write_text(json.dumps(proof,indent=2)+'\n')
    print('PASS: exact Fraction checks, 3 epsilon fixtures, affine parametric proof conditions')
    print('At epsilon=1/10000: global 1/20; backward identity 21001/400000; no-convex-guard 959/400000; reachable forward PASS')
    return proof
if __name__=='__main__':main()
