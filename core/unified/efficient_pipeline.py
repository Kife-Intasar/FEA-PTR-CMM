"""A SINGLE certificate state-machine
"""
from .exactbox import successor,residual
from .pipeline import verify_pair

def check_direct(source, rep, X0, X, U, guard):
    if guard is None:return {'status':'NO_GUARD'}
    cap=X0.within(guard) and guard.within(X)
    succ=successor(source,guard,U).intersect(X)
    closed=succ is None or succ.within(guard)
    resid=residual(source,rep,guard,U).within(rep.W)
    return {'status':'PASS' if cap and closed and resid else 'INCONCLUSIVE',
            'initial_capture':cap,'safe_successor_closure':closed,
            'representative_residual':resid}

def certify(source, rep, X0, X, U, H, candidate_guards, direct_guard=None,
            exact_lp_audit=False):
    guards=list(candidate_guards)
    if direct_guard is not None and direct_guard not in guards:guards.insert(0,direct_guard)
    for c in ([direct_guard] if direct_guard is not None else guards):
        d=check_direct(source,rep,X0,X,U,c)
        if d['status']=='PASS':
            return {'source':source.name,'representative':rep.name,'status':'CERTIFIED_DIRECT',
                    'finite_R1_T1':'LOGICALLY_ENTAILED','lp_facet_calls':0,
                    'initial_guard':d,'tail_B1_B2_B3':'ENTAILMENT'}
    finite=verify_pair(source,rep,X0,X,U,H,guards if guards else [None],run_lp=exact_lp_audit)
    if finite['status']=='CERTIFIED':
        empty=(finite.get('tail') or {}).get('status')=='CERTIFIED_EMPTY_TAIL'
        return {'source':source.name,'representative':rep.name,
                'status':'CERTIFIED_EMPTY_TAIL' if empty else 'CERTIFIED_TAIL',
                'finite_R1_T1':finite['finite']['status'],
                'lp_facet_calls':len(finite['lp']),
                'tail':finite.get('tail'),'evidence':'rational interval sufficient enclosure; LP optional numerical corroboration'}
    return {'source':source.name,'representative':rep.name,'status':'UNKNOWN',
            'lp_facet_calls':len(finite.get('lp',[])),
            'legacy_stage':finite['stage']}
