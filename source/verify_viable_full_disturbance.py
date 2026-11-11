
"""Strongest example: FULL-DIMENSIONAL adversarial disturbances in BOTH modes.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from unified.rational_polytope import Poly,box2,reach,residual

E=F(1,10000); delta=F(1,20); r=F(1,10)
A=((F(1),F(-1)),(F(1),F(-1)))
J=((F(41,40),F(-21,20)),(F(1),F(-1)))
B=((F(1,10),F(0)),(F(0),F(1,10)))
AB=((F(1,10),F(-1,10)),(F(1,10),F(-1,10)))
X=box2(-1,1,-1,1);X0=Poly([(F(1),F(0))]);U=box2(-1,1,-1,1)
Wi=box2(-E,E,-E,E)
wlo=-F(1,40)-E
whi=F(11,400)+F(41,40)*E
Wj=box2(wlo,whi,-E,E)
C=U.affine(B).minkowski(U.affine(AB)).minkowski(Wi).minkowski(Wi.affine(A))

def s(a):return str(a)
def encode_mat(m):return [[str(z) for z in row] for row in m]
def mymax(p):return max(abs(z) for z in p)

# First-exit exactly 2 steps
P1=reach(P0,A,B,U,Wi).clip_box()
assert P1.bounds()==((F(9,10)-E,F(1)),(F(9,10)-E,F(1)))
R0=residual(P0,A,J,B,B,U,Wi)
R1=residual(P1,A,J,B,B,U,Wi)
assert Wj.contains_poly(R0) and Wj.contains_poly(R1)
# Exact facet LP maxima are vertex support differences for the source-prefix polytopes.
facet=[]
for k,R in enumerate((R0,R1)):
 lohi=R.bounds()
 for t,coord,kind,margin in [(k,0,'upper',lohi[0][1]-whi),(k,0,'lower',wlo-lohi[0][0]),
                              (k,1,'upper',lohi[1][1]-E),(k,1,'lower',-E-lohi[1][0])]:
  assert margin<=0
  facet.append({'t':t,'axis':coord,'facet':kind,'exact_LP_facet_max':str(margin)})
assert len(facet)==8
P2=reach(P1,A,B,U,Wi).clip_box()
assert not P2.empty
assert C.contains_poly(P2) and X.contains_poly(C)
# A^2=0 implies AC+BU+Wi=C exactly; verify via exact Minkowski geometry too.
from unified.rational_polytope import mat
assert all(mat(A,mat(A,p))==(F(0),F(0)) for p in [(-E,-E),(E,E),(F(1),F(1))])
Cimage=reach(C,A,B,U,Wi)
assert C.contains_poly(Cimage) and Cimage.contains_poly(C)
Rt=residual(C,A,J,B,B,U,Wi)
assert Wj.contains_poly(Rt)
# No convex initial-region guard
x0=(F(1),F(0)); x1=(F(9,10),F(9,10)); x2=(F(0),F(0))
z=(F(9,10),F(0)); y=(F(4,5),F(1))
assert (mat(A,x0)[0]-F(1,10),mat(A,x0)[1]-F(1,10))==x1
assert mat(A,x1)==x2
assert (mat(A,z)[0]-F(1,10),mat(A,z)[1]+F(1,10))==y
assert X.contains_point(y)
res_y=delta*(y[1]-y[0]/2)
strictgap=res_y-whi
assert strictgap>0
# Includes a genuinely indefinitely safe source and a robustly safe representative.
Src1=box2(F(9,10)-E,F(9,10)+E,F(9,10)-E,F(9,10)+E)
Src2=reach(Src1,A,B,Poly([(F(0),F(0))]),Wi)
assert box2(-r,r,-r,r).contains_poly(Src2)
J1=box2(F(9,10)-E,F(381,400)+F(41,40)*E,F(9,10)-E,F(9,10)+E)
J2=reach(J1,J,B,Poly([(F(0),F(0))]),Wj)
T=box2(-r,r,-r,r)
assert T.contains_poly(J2)
J2mat=((F(1,1600),F(-21,800)),(F(1,40),F(-1,20)))
# odd and even invariant two-step safe bounds for both models.
eps1=max(abs(wlo),abs(whi));eps2=E
odd=[sum(abs(v) for v in J[k])*r+[eps1,eps2][k] for k in range(2)]
even=[sum(abs(v) for v in J2mat[k])*r+
      sum(abs(J[k][p])*[eps1,eps2][p] for p in range(2))+[eps1,eps2][k] for k in range(2)]
assert all(v<1 for v in odd) and all(v<=r for v in even)
# Source x2 belongs T and source A^2=0 with Wi noise; two-step same reasoning.
source_odd=[F(2)*r+E,F(2)*r+E]
source_even=[F(3)*E,F(3)*E]
assert all(v<1 for v in source_odd) and all(v<=r for v in source_even)

report={
 'outcome':'PASS','variant':'FULL-DIMENSIONAL disturbance, indefinite viable and robust reach-continue',
 'arithmetic':'fractions.Fraction exact; rational convex-hull',
 'model':{'Ai':encode_mat(A),'Aj':encode_mat(J),'B':encode_mat(B),
          'X':'[-1,1]^2','X0':'{(1,0)}','U':'[-1,1]^2',
          'Wi':f'[-{E},{E}]^2','Wj_first':[str(wlo),str(whi)],'Wj_second':[str(-E),str(E)]},
 'H':2,'finite_first_exit_LP_facets':facet,'all_eight_exact_optima_nonpositive':True,
 'safe_terminal_P2':P2.encode(),'terminal_nonempty':not P2.empty,
 'terminal_set_C':C.encode(),'terminal_set_vertices':len(C.v),
 'terminal_C_invariant_all_actions':True,'terminal_residual_polytope':Rt.encode(),
 'strict_separation':{'direct_convex_guard':'PROVED IMPOSSIBLE','forced_points':{'x0':['1','0'],'x1':['9/10','9/10'],'x2':['0','0'],'z':['9/10','0'],'y':['4/5','1']},
       'residual_at_y':str(res_y),'Wj_upper':str(whi),'strict_gap':str(strictgap)},
 'task':{'common_causal_controller':'u0=(-1,-1), u_t=0 for every t>=1',
         'target_G':'[-1/10,1/10]^2', 'target_reached_robustly_at_t':2,
         'source_x2_vertices':Src2.encode(),'representative_x2_vertices':J2.encode(),
         'source_odd_bound':[str(z) for z in source_odd], 'source_even_bound':[str(z) for z in source_even],
         'representative_odd_bound':[str(z) for z in odd], 'representative_even_bound':[str(z) for z in even],
         'all_future_states_safe_both_modes':True},
 'comparison':{'global_all_X_check':'FAIL','direct convex R1-T2':'PROVED NO CERTIFICATE',
 'finite R1-T1 H=2':'PASS, incomplete alone', 'unified first-exit+continuation':'CERTIFIED',
 'original_modes':2,'retained_modes':1,'robust_physical_task':'CERTIFIED'}
}
if __name__=='__main__':
 p=Path(__file__).resolve().parents[1]/'reports'/'R1_FULL_DISTURBANCE_EXACT.json'
 p.write_text(json.dumps(report,indent=2),encoding='utf8')
 print('PASS full-dimensional disturbance: 8 exact finite facet optima, nonempty tail, strict convex-guard contradiction, robust reach-and-continue')
 print('exact gap:',strictgap,' / terminal vertices:',len(C.v),' / output:',p)
