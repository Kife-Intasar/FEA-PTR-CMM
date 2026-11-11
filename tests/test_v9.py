import unittest
from fractions import Fraction as F
from global_mimo.exact_poly import Polytope,components
from global_mimo.history_bound import AffineMode,bound_at_history,combine_bounds
from global_mimo.scientific_run import basic_general_theory_tests,history_sensitive_test

class GeneralMIMOTheoremTests(unittest.TestCase):
    def test_geometry(self):self.assertEqual(basic_general_theory_tests()['general_3D_successor_components'],3)
    def test_history(self):self.assertEqual(history_sensitive_test()['history_based_general_lower_bound'],2)
    def test_poly_empty(self):
        with self.assertRaises(ValueError):Polytope.make([[1]],[F(-1)],[0],[1])
    def test_poly_point_touch(self):
        self.assertTrue(Polytope.box([0,0],[1,1]).intersects(Polytope.box([1,1],[2,2])))
    def test_poly_gap(self):
        self.assertFalse(Polytope.box([0,0],[1,1]).intersects(Polytope.box([F(1001,1000),0],[2,2])))
    def test_bad_posterior(self):
        W=Polytope.box(['-1/10'],['1/10'])
        m=AffineMode.make('z',[[0]],[[0]],W)
        with self.assertRaisesRegex(ValueError,'infeasible'):
            bound_at_history({'z':m},[(F(0),),(F(1),)],[(F(0),)],(F(0),))
    def test_mismatched_upper_rejected(self):
        with self.assertRaisesRegex(ValueError,'Contradictory'):
            combine_bounds([{'global_mode_lower_bound':3}],1)
    def test_unknown_no_upper(self):
        r=combine_bounds([{'global_mode_lower_bound':1}])
        self.assertIsNone(r['globally_minimal'])
    def test_inconclusive_gap(self):
        r=combine_bounds([{'global_mode_lower_bound':2}],3)
        self.assertEqual(r['status'],'GAP_OPEN')
    def test_disconnected_budget(self):
        W=Polytope.box([0],[1])
        mods={'a':AffineMode.make('a',[[0]],[[1]],W),
              'b':AffineMode.make('b',[[0]],[[1]],W.shifted([3]))}
        r=bound_at_history(mods,[(F(0),)],[],(F(0),),2)
        self.assertEqual(r['global_mode_lower_bound'],1)
    def test_mode_history_not_switched(self):
        self.assertEqual(history_sensitive_test()['separate_posterior_pruning_test'],['plus'])

if __name__=='__main__':unittest.main()

class RepresentationBoundaryTests(unittest.TestCase):
    def test_memory_lift(self):
        from global_mimo.memory_lift_counterexample import audit
        self.assertEqual(audit()['result'],'PASS')

class NonlinearMIMOCertificationTests(unittest.TestCase):
    def test_nonlinear_coupled_global_minimum(self):
        from global_mimo.nonlinear_mimo_fixture import run
        self.assertEqual(run()['global_minimum_same_physical_state_connected_images'],3)

class CostGapGeneralisationTests(unittest.TestCase):
    def test_connected_outer_sharp(self):
        from global_mimo.connected_cost_gap import sharp_connected_outer_gap
        d=sharp_connected_outer_gap(F(7,100),F(13,100),F(1,200),F(3,200))
        self.assertEqual(d['sharp_minimum_worst_case_cost_of_ANY_CONNECTED_OUTER_APPROXIMATION'],'7/100')
    def test_no_gap_if_touch(self):
        from global_mimo.connected_cost_gap import sharp_connected_outer_gap
        with self.assertRaises(ValueError):sharp_connected_outer_gap(0,1,0,1)
