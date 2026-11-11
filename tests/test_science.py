import unittest
from fractions import Fraction as Q
from unified.exactbox import Box,Mode,finite_outer,continuation,residual,successor
from unified.fixtures import fourteen,bridge,strict
from unified.pipeline import verify_pair,cover,verify_exact_counterexample

class ScientificChecks(unittest.TestCase):
    def test_exact_negative_witness(self):
        m,p,x0,X,U,C=fourteen()
        source=m[p[0][0]];wrong=m['R-1' if source.name.startswith('I+') else 'R+1']
        z=verify_exact_counterexample(source,wrong,(Q(0),Q(0)),[(Q(1,10),Q(0))],[(Q(0),Q(0))],X,U)
        self.assertEqual(z['status'],'REJECT_EXACT')
    def test_fourteen(self):
        mods,pairs,x0,X,U,C=fourteen()
        good=[]
        for a,b in pairs:
            z=verify_pair(mods[a],mods[b],x0,X,U,4,[C],run_lp=False)
            self.assertEqual(z['status'],'CERTIFIED')
            good.append((a,b))
        self.assertEqual(len(cover(list(mods),good)),2)
    def test_bridge_nonempty_terminal(self):
        i,j,X0,X,U,C=bridge()
        z=verify_pair(i,j,X0,X,U,3,[C],run_lp=False)
        self.assertEqual(z['status'],'CERTIFIED')
        self.assertIsNotNone(z['finite']['safe_terminal'])
        self.assertEqual(z['tail']['status'],'CERTIFIED')
    def test_strict_finite_exit(self):
        i,j,X0,X,U=strict()
        z=verify_pair(i,j,X0,X,U,3,[None],run_lp=False)
        self.assertIsNone(z['finite']['safe_terminal'])
        self.assertEqual(z['tail']['status'],'CERTIFIED_EMPTY_TAIL')
    def test_inconclusive_is_not_impossible(self):
        i,j,X0,X,U,C=bridge()
        z=verify_pair(i,j,X0,X,U,3,[],run_lp=False)
        self.assertEqual(z['status'],'UNKNOWN')
    def test_first_exit_includes_unsafe_successor(self):
        i=Mode.make('a',[[0]],[['0']],Box.singleton([3]))
        j=Mode.make('b',[[0]],[['0']],Box.singleton([0]))
        z=finite_outer(i,j,Box.singleton([0]),Box.symmetric([1]),Box.singleton([0]),1)
        self.assertEqual(z['status'],'UNKNOWN_OUTER')
    def test_direct_implies_finite(self):
        m,p,x0,X,U,C=fourteen()
        i,j=p[0]
        self.assertTrue(successor(m[i],C,U).intersect(X).within(C))
        self.assertTrue(residual(m[i],m[j],C,U).within(m[j].W))
        self.assertEqual(finite_outer(m[i],m[j],x0,X,U,4)['status'],'PASS')

if __name__=='__main__':unittest.main()
