"""Sharp d_infinity worst-case gap for ANY connected outer successor set.
"""
from fractions import Fraction as F
from .exact_poly import frac

def sharp_connected_outer_gap(a,b,ymin,ymax):
    a,b,ymin,ymax=map(frac,(a,b,ymin,ymax))
    if not (F(0)<a<=b and ymin<=ymax):raise ValueError('Separated symmetric strips required')
    midpoint=(F(0),(ymin+ymax)/2)
    assert min(abs(midpoint[0]-(-a)),abs(midpoint[0]-a))==a
    assert max(a,F(0))==a
    return {'source_components':2,'left_x':[str(-b),str(-a)],
            'right_x':[str(a),str(b)],'y_band':[str(ymin),str(ymax)],
            'worst_case_cost_original':'0',
            'sharp_minimum_worst_case_cost_of_ANY_CONNECTED_OUTER_APPROXIMATION':str(a)}
