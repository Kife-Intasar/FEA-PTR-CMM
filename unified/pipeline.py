
from .exactbox import finite_outer, continuation
from .facetlp import solve
from itertools import combinations

def verify_pair(i,j,X0,X,U,H,candidate_guards,run_lp=True):
    finite=finite_outer(i,j,X0,X,U,H)
    numerical=solve(i,j,X0,X,U,H) if run_lp else []
    if finite['status']=='UNKNOWN_OUTER':
        return {'source':i.name,'representative':j.name,'H':H,'status':'UNKNOWN','stage':'finite_outer',
                'finite':serial(finite),'lp':numerical}
    for G in candidate_guards:
        tail=continuation(i,j,finite['safe_terminal'],X,U,G)
        if tail['status'].startswith('CERTIFIED'):
            return {'source':i.name,'representative':j.name,'H':H,'status':'CERTIFIED','stage':'continuation',
                    'finite':serial(finite),'lp':numerical,'tail':tail,
                    'reason':'sound exact Fraction interval certificate; not rational LP optimum'}
    return {'source':i.name,'representative':j.name,'H':H,'status':'UNKNOWN',
            'stage':'continuation','finite':serial(finite),'lp':numerical}

def serial(z):
    v=dict(z)
    if v.get('safe_terminal') is not None:v['safe_terminal']=v['safe_terminal'].export()
    return v

def cover(names,arcs):
    for q in range(1,len(names)+1):
        for C in combinations(names,q):
            kept=set(C)
            if all(i in kept or any((i,j) in arcs for j in kept) for i in names):return list(C)
    return list(names)

def verify_exact_counterexample(source,representative,x0,controls,disturbances,X,U):
    from .exactbox import Box,linear,add
    def singleton(a):return Box.singleton(a)
    if not singleton(x0).within(X):return {'status':'INVALID_WITNESS'}
    x=singleton(x0)
    if len(controls)!=len(disturbances) or not controls:return {'status':'INVALID_WITNESS'}
    for t,(u,w) in enumerate(zip(controls,disturbances)):
        ub=singleton(u);wb=singleton(w)
        if not ub.within(U) or not wb.within(source.W):return {'status':'INVALID_WITNESS'}
        y=add(linear(source.A,x),linear(source.B,ub),wb)
        inferred=add(y,linear(tuple(tuple(-p for p in row) for row in representative.A),x),
                        linear(tuple(tuple(-p for p in row) for row in representative.B),ub))
        if not inferred.within(representative.W):
            return {'status':'REJECT_EXACT','first_failing_transition':t,'target_residual':inferred.export(),
                    'unsafe_successor_allowed':True}
        x=y
        if t+1<len(controls) and not x.within(X):return {'status':'INVALID_WITNESS'}
    return {'status':'NO_VIOLATION_IN_WITNESS'}
