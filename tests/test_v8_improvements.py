import unittest
from fractions import Fraction as F
from unified.fixtures import fourteen,bridge
from unified.global_minimum import exact_audit
from unified.efficient_pipeline import certify
from unified.pipeline import verify_exact_counterexample
class V8Tests(unittest.TestCase):
    def test_global_min(self):
        r=exact_audit();self.assertEqual(r['global_minimum_universal_worst_case_physical_costs'],2)
        self.assertEqual(r['candidate_original_subsets_exhausted'],16383)
        self.assertEqual(r['single_prescribed_X_safety_truth_minimum'],1)
        self.assertEqual(r['exact_one_step_conservative_one_mode_cost_gap_at_least'],'7/100')
    def test_guard_first_all_12(self):
        mods,pairs,x0,X,U,C=fourteen()
        for a,b in pairs:
            s=certify(mods[a],mods[b],x0,X,U,4,[C],direct_guard=C)
            self.assertEqual(s['status'],'CERTIFIED_DIRECT')
            self.assertEqual(s['lp_facet_calls'],0)
            self.assertEqual(s['finite_R1_T1'],'LOGICALLY_ENTAILED')
    def test_terminal_only_when_initial_guard_not_available(self):
        i,j,x0,X,U,C=bridge()
        s=certify(i,j,x0,X,U,3,[C],direct_guard=None)
        self.assertIn(s['status'],('CERTIFIED_TAIL','CERTIFIED_DIRECT'))
    def test_empty_guard_unknown_not_impossible(self):
        i,j,x0,X,U,C=bridge()
        s=certify(i,j,x0,X,U,3,[],direct_guard=None)
        self.assertEqual(s['status'],'UNKNOWN')
    def test_negative_wrong_group_exact_witness(self):
        modes,pairs,x0,X,U,C=fourteen()
        s=modes[pairs[0][0]]
        wrong=modes['R-1' if s.name.startswith('I+') else 'R+1']
        out=verify_exact_counterexample(s,wrong,(F(0),F(0)),[(F(1,10),F(0))],[(F(0),F(0))],X,U)
        self.assertEqual(out['status'],'REJECT_EXACT')
if __name__=='__main__':unittest.main()
