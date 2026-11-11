"""Exact global minimum on the fourteen-mode fixture, with explicit quantifier scope.
"""
from fractions import Fraction as F
from itertools import combinations
from .fixtures import fourteen
from .exactbox import Box, linear, add, successor

def one_step_rectangles():
    modes, pairs, x0, X, U, C = fourteen()
    assert x0.lo == (F(0),F(0)) and x0.hi == (F(0),F(0))
    u=Box.singleton((F(1,10),F(0)))
    assert u.within(U)
    out={}
    for name,mode in modes.items():
        z=add(linear(mode.A,x0),linear(mode.B,u),mode.W)
        assert z.within(X)  
        out[name]=z
    return out

def intersect(a,b):
    return all(max(a.lo[k],b.lo[k])<=min(a.hi[k],b.hi[k]) for k in range(2))

def components(boxes):
    keys=list(boxes)
    adj={k:set() for k in keys}
    for a,b in combinations(keys,2):
        if intersect(boxes[a],boxes[b]):
            adj[a].add(b);adj[b].add(a)
    unseen=set(keys);parts=[]
    while unseen:
        cur={next(iter(unseen))};seen=set(cur)
        while cur:
            nxt=set().union(*(adj[s] for s in cur))-seen
            seen.update(nxt);cur=nxt
        unseen-=seen;parts.append(seen)
    return parts

def point_in_box(p,b):
    return all(b.lo[k]<=p[k]<=b.hi[k] for k in (0,1))

def exact_audit():
    rect=one_step_rectangles();comps=components(rect)
    modes,_,X0,X,U,C=fourteen()
    assert X0.within(C) and C.within(X)
    assert all(successor(m,C,U).within(C) for m in modes.values())
    assert len(rect)==14 and len(comps)==2
    left=rect['R-1'];right=rect['R+1']
    assert left.lo==(F(-13,100),F(1,200)) and left.hi==(F(-7,100),F(3,200))
    assert right.lo==(F(7,100),F(1,200)) and right.hi==(F(13,100),F(3,200))
    for key,b in rect.items():
        target=right if key.startswith(('R+','I+')) else left
        assert b.within(target)
    assert not intersect(left,right)
    edges=sum(intersect(a,b) for a,b in combinations(rect.values(),2))
    assert edges==42
    qleft=(F(-13,100),F(1,100));qright=(F(13,100),F(1,100))
    assert {s for s,b in rect.items() if point_in_box(qleft,b)}=={'R-1'}
    assert {s for s,b in rect.items() if point_in_box(qright,b)}=={'R+1'}
    names=tuple(rect);winning=[];tested=0
    for mask in range(1,1<<len(names)):
        tested+=1
        ss={names[k] for k in range(len(names)) if (mask>>k)&1}
        endpoints=any(point_in_box(qleft,rect[k]) for k in ss) and any(point_in_box(qright,rect[k]) for k in ss)
        assert endpoints==({'R+1','R-1'}<=ss)
        if endpoints and len(ss)==2:winning.append(sorted(ss))
    assert tested==16383 and winning==[['R+1','R-1']]
    m=(F(0),F(1,100))
    assert not (point_in_box(m,left) or point_in_box(m,right))
    lower=min(abs(m[0]-left.hi[0]),abs(m[0]-right.lo[0]))
    assert lower==F(7,100)
    return {
      'audit':'PASS','arithmetic':'Python fractions; exact',
      'original_modes':14,'exact_original_successor_components':2,
      'within_group_intersecting_pairs':42,'cross_group_disjoint_pairs':49,
      'one_step_probe_x':'(0,0)','one_step_probe_u':'(1/10,0)',
      'left_box':left.export(),'right_box':right.export(),
      'candidate_original_subsets_exhausted':tested,
      'optimal_two_subset':winning[0],
      'minimum_over_all_original_mode_subsets_exact_one_step':2,
      'global_lower_bound_over_any_single_affine_convex_additive_mode':2,
      'global_upper_bound_from_frozen_V7_trace_certificate':2,
      'global_minimum_same_state_convex_additive_affine_universal_physical_traces':2,
      'global_minimum_universal_worst_case_physical_costs':2,
      'all_modes_or_arbitrary_nonconvex_W':'NOT PROVED',
      'single_prescribed_X_safety_truth_minimum':1,
      'exact_one_step_conservative_one_mode_cost_gap_at_least':str(lower),
      'explicit_cost':'d_infty(x_1, union_of_original_probe_successors)'
    }
