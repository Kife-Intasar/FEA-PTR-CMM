import unittest,json
from pathlib import Path
from fractions import Fraction as Q
from unified.viable_unified import pipeline,data
from unified.rational_polytope import box2,Poly,reach,residual

class ViableStrict(unittest.TestCase):
 def test_exact_full_dim_certification(self):
  a=pipeline(2,True)
  self.assertEqual(a['status'],'CERTIFIED')
  self.assertTrue(a['terminal_safe_nonempty'])
  self.assertEqual(len(a['finite_stages']),2)
  self.assertTrue(all(x['representative_facet_pass'] for x in a['finite_stages']))
  self.assertTrue(a['tail']['all_action_invariant'])
  self.assertTrue(a['tail']['residual_included'])
 def test_full_dim_disturbances(self):
  _,_,_,_,_,_,Wi,Wj,_=data(True)
  self.assertGreater(Wi.bounds()[0][1]-Wi.bounds()[0][0],0)
  self.assertGreater(Wi.bounds()[1][1]-Wi.bounds()[1][0],0)
  self.assertGreater(Wj.bounds()[0][1]-Wj.bounds()[0][0],0)
  self.assertGreater(Wj.bounds()[1][1]-Wj.bounds()[1][0],0)
 def test_nonempty_terminal_exact_clipping(self):
  A,J,B,X,P,U,Wi,Wj,C=data(True)
  self.assertFalse(reach(reach(P,A,B,U,Wi).clip_box(),A,B,U,Wi).clip_box().empty)
 def test_non_convex_failure_point(self):
  A,J,B,X,P,U,Wi,Wj,C=data(True)
  y=(Q(4,5),Q(1))
  self.assertTrue(X.contains_point(y))
  D=((Q(-1,2),Q(1)),(Q(0),Q(0)))
  wr=Q(1,20)*(y[1]-y[0]/2)
  self.assertGreater(wr,Wj.bounds()[0][1])
 def test_original_r1_t2_not_wholly_replaced(self):
  A,J,B,X,P,U,Wi,Wj,C=data(True)
  self.assertTrue(X.contains_poly(C))
 def test_degenerate_zero_disturbance_special_case(self):
  self.assertEqual(pipeline(2,False)['status'],'CERTIFIED')

if __name__=='__main__': unittest.main()
