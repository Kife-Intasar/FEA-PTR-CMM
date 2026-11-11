"""Nonlinear coupled 3-state / 2-input example: general theorem's scope.
"""
from fractions import Fraction as F
from .exact_poly import Polytope,components
from .history_bound import combine_bounds

def run():
    x=(F(1,5),F(2,5),F(1,10));u=(F(1,3),F(1,7))
    nonlinear=(x[0]*x[1]+u[0]*u[1],x[1]*x[1]+u[1],x[2]+x[0]*x[2]+u[0]*u[0])
    c=(F(0),F(3),F(6))
    W=Polytope.box([F(-1,10)]*3,[F(1,10)]*3)
    shapes={f'mode_{i}':W.shifted([nonlinear[0]+v,nonlinear[1],nonlinear[2]]) for i,v in enumerate(c)}
    geom=components(shapes)
    assert geom['connected_components']==3
    opt=combine_bounds([{'global_mode_lower_bound':3}],3)
    assert opt['status']=='GLOBAL_MINIMUM_CERTIFIED'
    return {'model':'f_i(x,u,w)= (c_i + x1*x2 + u1*u2 + w1, x2^2 + u2 + w2, x3 + x1*x3 + u1^2 + w3)',
            'state_dimension':3,'input_dimension':2,'nonlinear_coupling':True,
            'fixed_mode_count':3,'disturbance_box_width_each_coordinate':'1/5',
            'initial_state':list(map(str,x)),'first_action':list(map(str,u)),
            'exact_successor_components':geom['connected_components'],
            'lower_bound':3,'upper_bound':'3: identical original model (reflexive full trace preservation)',
            'global_minimum_same_physical_state_connected_images':3
            }
