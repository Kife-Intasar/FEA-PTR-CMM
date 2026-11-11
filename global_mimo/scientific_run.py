
from pathlib import Path
from fractions import Fraction as F
import json,sys,tempfile,zipfile,subprocess,os
from .exact_poly import Polytope, components
from .history_bound import AffineMode,bound_at_history,combine_bounds
from .nonlinear_mimo_fixture import run as nonlinear_audit
from .memory_lift_counterexample import audit as memory_audit
from .connected_cost_gap import sharp_connected_outer_gap


def basic_general_theory_tests():
    diamond=Polytope.make([[1,1],[-1,1],[1,-1],[-1,-1]],[1,1,1,1],[-1,-1],[1,1])
    right=Polytope.box([1,0],[1,0])
    away=Polytope.box([F(101,100),0],[2,1])
    assert diamond.intersects(right)  # touches at exactly (1,0)
    assert not diamond.intersects(away)
    c=components({'oblique':diamond,'boundary_tangent':right,'distant':away})
    assert c['connected_components']==2 and c['intersecting_pairs']==1
    tilted=Polytope.make([[1,1,1]],[F(3,2)],[-1,-1,-1],[1,1,1])
    near=Polytope.box([1,0,0],[F(3,2),F(1,2),F(1,2)])
    far=Polytope.box([3,3,3],[4,4,4])
    g=components({'tilted':tilted,'near':near,'far':far})
    assert g['connected_components']==2
    multi=components({'a':Polytope.box([0,0,0],[1,1,1]),
                      'b':Polytope.box([3,0,0],[4,1,1]),
                      'c':Polytope.box([6,0,0],[7,1,1])})
    assert multi['connected_components']==3
    three=combine_bounds([{'global_mode_lower_bound':3}],validated_upper=3)
    assert three['globally_minimal'] and three['lower_bound']==3
    unresolved=combine_bounds([{'global_mode_lower_bound':2}],validated_upper=3)
    assert unresolved['status']=='GAP' and not unresolved['globally_minimal']
    one=combine_bounds([{'global_mode_lower_bound':1}])
    assert one['status']=='LOWER_BOUND'
    budget=combine_bounds([{'global_mode_lower_bound':(3+2)//3}],validated_upper=1)
    assert budget['globally_minimal']
    return {'general_3D_successor_components':3,
            'general_3D_global_minimum_by_original_identity_upper_bound':3,
            'disconnected_per_mode_successor_budget_test':1}


def history_sensitive_test():
    W=Polytope.box(['-1/10','-1/10'],['1/10','1/10'])
    Z=[['0','0'],['0','0']]
    mi={
      'plus':AffineMode.make('plus',Z,[['1','0'],['0','1']],W),
      'minus':AffineMode.make('minus',Z,[['-1','0'],['0','1']],W)
    }
    initial=bound_at_history(mi,[(F(0),F(0))],[],(F(0),F(0)))
    later=bound_at_history(mi,[(F(0),F(0)),(F(0),F(0))],[(F(0),F(0))],(F(1),F(0)))
    assert initial['global_mode_lower_bound']==1
    assert later['global_mode_lower_bound']==2
    assert later['posterior_modes']==['minus','plus']
    prune=bound_at_history(mi,[(F(0),F(0)),(F(1),F(0))],[(F(1),F(0))],(F(0),F(0)))
    assert prune['posterior_modes']==['plus'] and prune['global_mode_lower_bound']==1
    return {'two_mode_MIMO_input_dim':2,'physical_state_dim':2,
            'initial_selected_probe_components':1,
            'history_conditioned_components':2,
            'history_based_general_lower_bound':2,
            'fixed_mode_posterior_test':'PASS',
            'history_posterior':['minus','plus'],
            'separate_posterior_pruning_test':['plus']}


def original_14_and_matching(V_root:Path):
    core=V_root/'core'
    if not core.is_dir():raise FileNotFoundError(core)
    sys.path.insert(0,str(core))
    try:
        from unified.fixtures import fourteen
        from unified.efficient_pipeline import certify
        from unified.global_minimum import exact_audit
        old=exact_audit()
        modes,pairs,X0,X,U,C=fourteen()
        assert len(modes)==14 and len(pairs)==12
        x=X0.lo; u=(F(1,10),F(0))
        transformed={name:AffineMode.make(name,m.A,m.B,
           Polytope.box(m.W.lo,m.W.hi)) for name,m in modes.items()}
        probe=bound_at_history(transformed,[x],[],u)
        assert probe['successor_connected_components']==2
        assert probe['geometry']['intersecting_pairs']==42
        assert probe['geometry']['disjoint_pairs']==49
        results=[certify(modes[i],modes[j],X0,X,U,4,[C],direct_guard=C) for i,j in pairs]
        assert all(z['status']=='CERTIFIED_DIRECT' and z['lp_facet_calls']==0 for z in results)
        retained=sorted(set(j for _,j in pairs))
        assert retained==['R+1','R-1']
        upper=len(retained)
        match=combine_bounds([probe],validated_upper=upper)
        assert match['status']=='GLOBAL_MINIMUM_CERTIFIED' and match['lower_bound']==2
        assert old['global_minimum_same_state_convex_additive_affine_universal_physical_traces']==2
        assert old['single_prescribed_X_safety_truth_minimum']==1
        return {'status':match['status'],
                'original_modes':14,'certified_retained_modes':retained,
                'minimum':2, 'topological_lower_bound':2,'independent_U_T1_upper_bound':upper,
                'matched_certificate_pairs':len(results), 'original_subset_enumeration':old['candidate_original_subsets_exhausted'],
                'probe':{'x':list(map(str,x)),'u':list(map(str,u))},
                'geometry':probe['geometry'],
                'one_prescribed_X_safety_task_minimum':1
               }
    finally:
        sys.path.remove(str(core))


def run(V_root):
    return {
       'general_exact_geometry':basic_general_theory_tests(),
       'coupled_nonlinear_3D_2input_global_minimum':nonlinear_audit(),
       'unrestricted_memory_lift_counterexample':memory_audit(),
       'sharp_connected_outer_cost_gap':sharp_connected_outer_gap('7/100','13/100','1/200','3/200'),
       'history_conditioning':history_sensitive_test(),
       'original_14_mode_global_mimo_certificate':original_14_and_matching(V_root)}


if __name__=='__main__':
    V=Path(sys.argv[1]) if len(sys.argv)>1 else None
    if V is None:raise SystemExit('Usage: python -m global_mimo.scientific_run')
    print(json.dumps(run(V),indent=2))
