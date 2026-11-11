
from fractions import Fraction as Q
from itertools import combinations
from .fixtures import fourteen

def geometry():
    mods,pairs,X0,X,U,C=fourteen();u=Q(1,10);r=[]
    for name,m in mods.items():
        center=[m.B[k][0]*u for k in range(2)]
        r.append((name,1 if name.startswith(('R+','I+')) else -1,
                 [center[k]+m.W.lo[k] for k in range(2)],
                 [center[k]+m.W.hi[k] for k in range(2)]))
    def intersection(a,b):return all(max(a[2][k],b[2][k])<=min(a[3][k],b[3][k]) for k in (0,1))
    within=cross=0
    for a,b in combinations(r,2):
        if a[1]==b[1]:assert intersection(a,b);within+=1
        else:assert not intersection(a,b);cross+=1
    assert (within,cross)==(42,49)
    negative=[a for a in r if a[1]==-1];positive=[a for a in r if a[1]==1]
    assert min(a[2][0] for a in positive)==Q(7,100)
    assert max(a[3][0] for a in negative)==-Q(7,100)
    return {'within_cluster_intersections':within,'cross_cluster_disjoint':cross,
            'connected_components':2,'minimum_convex_affine_exact_modes':2,
            'best_one_step_convex_Hausdorff_infty':'7/100',
            'warning':'fixed admissible probe and exact convex-additive representation class only'}

def terminal_decisions():
    universe=('a','b','c')
    subsets=[set(z) for k in range(4) for z in combinations(universe,k)]
    count=0
    for r in subsets:
        for omit in subsets:
            assert ((r & omit)==r)==r.issubset(omit);count+=1
    return {'checked_finite_set_intersection_instances':count,
            'conclusion':'D_original = D_retained iff all omitted original-consistent modes allow every retained-robust decision',
            'interpretation':'only verification of established set algebra, not newly discovered theorem'}

def absorption():
    lam=Q(3,5); c=Q(14,125);R=Q(13,10);r=Q(7,10)
    value=lambda k:lam**k*R+c*sum((lam**s for s in range(k)),Q(0))
    assert value(1)==Q(223,250) and value(2)==Q(809,1250)
    assert value(1)>r and value(2)<r
    assert lam*r+c<r
    return {'all_action_upper_bounds':[str(value(k)) for k in range(3)],
            'absorbing_radius':str(r),'certified_horizon':2,
            'restriction':'uniformly contracting special subclass; no general MIMO finite determination'}
