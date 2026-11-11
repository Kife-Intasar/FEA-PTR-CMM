"""Fresh method-stage timing, same candidates, identical rational fixture.
"""
from time import perf_counter
from statistics import median
from .fixtures import fourteen
from .exactbox import residual,successor
from .pipeline import verify_pair

def benchmark(repetitions=9):
    modes,pairs,x0,X,U,C=fourteen()
    direct=[];unified=[];interval=[];globaltest=[]
    for _ in range(repetitions):
        ts=perf_counter();nr=0
        for a,b in pairs:
            i=modes[a];j=modes[b]
            good=x0.within(C) and C.within(X) and successor(i,C,U).intersect(X).within(C) and residual(i,j,C,U).within(j.W)
            nr+=bool(good)
        direct.append(1000*(perf_counter()-ts));assert nr==12
        ts=perf_counter();nr=0
        for a,b in pairs:
            good=residual(modes[a],modes[b],X,U).within(modes[b].W)
            nr+=bool(good)
        globaltest.append(1000*(perf_counter()-ts));assert nr==0
        ts=perf_counter();nr=0
        for a,b in pairs:
            q=verify_pair(modes[a],modes[b],x0,X,U,4,[C],run_lp=False)
            nr+=(q['status']=='CERTIFIED')
        interval.append(1000*(perf_counter()-ts));assert nr==12
        ts=perf_counter();nr=0
        for a,b in pairs:
            q=verify_pair(modes[a],modes[b],x0,X,U,4,[C],run_lp=True)
            nr+=(q['status']=='CERTIFIED')
        unified.append(1000*(perf_counter()-ts));assert nr==12
    return {'stage':'12 fixed candidate source-representative pairs; supplied identical guard',
            'units':'milliseconds; measured fresh on present execution host; not general timing claim',
            'N_repetitions':repetitions,'global_check_median_ms':median(globaltest),
            'direct_guard_median_ms':median(direct),'unified_fraction_outer_plus_tail_median_ms':median(interval),
            'unified_finite_LP_plus_tail_median_ms':median(unified),
            'unified_fraction_vs_direct_ratio':median(interval)/median(direct),
            'unified_with_LP_vs_direct_ratio':median(unified)/median(direct),
            'all_12_global_fail':True, 'direct_guard_12_certified':True, 'unified_12_certified':True
          }
