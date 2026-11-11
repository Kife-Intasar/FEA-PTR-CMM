"""Fixed-mode Bayesian-free posterior and topological minimum bound.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import ceil
from .exact_poly import Polytope, frac, dot, components


def mat_vec(M,x):
    return tuple(dot(row,x) for row in M)


@dataclass(frozen=True)
class AffineMode:
    name:str
    A:tuple
    B:tuple
    W:Polytope
    @classmethod
    def make(cls,name,A,B,W):
        Aa=tuple(tuple(map(frac,r)) for r in A)
        Bb=tuple(tuple(map(frac,r)) for r in B)
        n=W.dim
        if len(Aa)!=n or any(len(r)!=n for r in Aa) or len(Bb)!=n:raise ValueError('state dimensions')
        if len({len(r) for r in Bb})!=1 or len(Bb[0])==0:raise ValueError('input dimensions')
        return cls(name,Aa,Bb,W)
    def displacement(self,x,u):
        ax=mat_vec(self.A,x);bu=mat_vec(self.B,u)
        return tuple(a+b for a,b in zip(ax,bu))
    def contains_step(self,x,u,y):
        z=self.displacement(x,u)
        return self.W.contains(tuple(frac(y[k])-z[k] for k in range(self.W.dim)))
    def successor(self,x,u):
        return self.W.shifted(self.displacement(x,u))


def posterior(modes,states,controls):
    if len(states)!=len(controls)+1 or not states:raise ValueError('History length')
    return {name:mode for name,mode in modes.items() if all(
        mode.contains_step(states[t],controls[t],states[t+1]) for t in range(len(controls)))}


def bound_at_history(modes,states,controls,next_u,max_components_per_replacement=1):
    post=posterior(modes,states,controls)
    if not post:raise ValueError('Provided physical history is infeasible for every fixed mode')
    shapes={name:mode.successor(states[-1],next_u) for name,mode in post.items()}
    info=components(shapes)
    c=int(max_components_per_replacement)
    if c<=0:raise ValueError('Every representative needs a positive component budget')
    return {'posterior_modes':sorted(post), 'history_steps':len(controls),
            'successor_connected_components':info['connected_components'],
            'candidate_mode_component_budget':c,
            'global_mode_lower_bound':ceil(info['connected_components']/c),
            'geometry':info}


def combine_bounds(results,validated_upper=None):
    if not results:raise ValueError('Supply at least one EXACT reachable-history certificate')
    lb=max(v['global_mode_lower_bound'] for v in results)
    out={'lower_bound':lb,'upper_bound':None,'status':'LOWER_BOUND','globally_minimal':None}
    if validated_upper is not None:
        if type(validated_upper)!=int or validated_upper<1:raise ValueError('Upper bound must be a verified positive integer')
        if validated_upper<lb:raise ValueError('Contradictory upper/lower certificates!')
        out.update(upper_bound=validated_upper,
                   status='GLOBAL_MINIMUM_CERTIFIED' if lb==validated_upper else 'GAP',
                   globally_minimal=(lb==validated_upper))
    return out
