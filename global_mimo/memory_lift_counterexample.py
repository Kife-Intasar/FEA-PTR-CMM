"""Executable finite-horizon sanity example for the unrestricted representation warning.
"""
from fractions import Fraction as F
from itertools import product

def actions():return ((F(0),F(0)),(F(1),F(0)))
def disturbances():return ((F(0),F(0)),(F(1,10),F(0)))
def step(mode,x,u,w):
    sign=1 if mode=='plus' else -1
    return (x[0]+sign*u[0]+w[0],x[1]+u[1]+w[1])

def original(T):
    ans=set()
    for mode in ('plus','minus'):
        for us in product(actions(), repeat=T):
            for ws in product(disturbances(), repeat=T):
                x=(F(0),F(0));trace=[x]
                for u,w in zip(us,ws):
                    x=step(mode,x,u,w);trace.append(x)
                ans.add((us,tuple(trace)))
    return ans

def augmented_one_rule(T):
    ans=set()
    for theta in ('plus','minus'):
        for us in product(actions(), repeat=T):
            for ws in product(disturbances(), repeat=T):
                z=((F(0),F(0)),theta);trace=[z[0]]
                for u,w in zip(us,ws):
                    z=(step(z[1],z[0],u,w),z[1]);trace.append(z[0])
                ans.add((us,tuple(trace)))
    return ans

def audit():
    counts=[]
    for T in (1,2,3):
        a=original(T);b=augmented_one_rule(T)
        assert a==b;counts.append({'T':T,'distinct_physical_input_trajectories':len(a)})
    return {'result':'PASS','original_fixed_modes':2,'explicit_mode_labels_in_augmented_governing_rule':1,
            'finite_horizon_checks':counts}
