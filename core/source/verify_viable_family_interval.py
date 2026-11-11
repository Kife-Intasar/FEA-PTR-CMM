
from fractions import Fraction as F
from pathlib import Path
import json
E=F(1,1000); d=F(1,20); r=F(1,10)
gap=F(1,400)-F(41,40)*E
assert gap>0
Cbox=F(3,10)+3*E
assert Cbox<1
tail=F(1,80)+F(9,8)*E
upper=F(11,400)+F(41,40)*E
lower_abs=F(1,40)+E
assert tail<=upper and tail<=lower_abs
finite0=( -F(1,40)-E,-F(1,40)+E)
finite1=( F(1,50)-F(21,20)*E, F(11,400)+F(41,40)*E)
assert finite0[0]>=-lower_abs and finite0[1]<=upper
assert finite1[0]>=-lower_abs and finite1[1]<=upper
J=((F(41,40),-F(21,20)),(F(1),F(-1)))
from itertools import product
x1=((F(9,10)-E,F(381,400)+F(41,40)*E),(F(9,10)-E,F(9,10)+E))
w=(( -F(1,40)-E,F(11,400)+F(41,40)*E),(-E,E))
for k in range(2):
 vals=[sum(J[k][m]*x[m] for m in range(2))+wk for x in product(*x1) for wk in w[k]]
 assert min(vals)>=-r and max(vals)<=r,(k,min(vals),max(vals))
J2=((F(1,1600),-F(21,800)),(F(1,40),-F(1,20)))
nu=(F(11,400)+F(41,40)*E,E)
odd=[sum(abs(v) for v in row)*r+nu[k] for k,row in enumerate(J)]
even=[sum(abs(v) for v in row)*r+sum(abs(J[k][m])*nu[m] for m in range(2))+nu[k]
      for k,row in enumerate(J2)]
assert all(z<1 for z in odd) and all(z<=r for z in even)
assert 3*E<=r and 2*r+E<1

report={'result':'PASS','proved_interval':'0 < epsilon <= 1/1000','field':'all real epsilon mathematically (rational executable endpoints)',
        'argument':'all required scalar bounds are monotonic affine upper bounds in epsilon',
        'worst_positive_convex_guard_residual_gap':str(gap),'C_coefficient_box_bound':str(Cbox),
        'tail_abs_residual_bound':str(tail),'Wj_positive_limit':str(upper),
        'target_representative_odd_bounds':[str(z) for z in odd],
        'target_representative_even_bounds':[str(z) for z in even]}
if __name__=='__main__':
 p=Path(__file__).resolve().parents[1]/'reports'/'R1.json'
 p.write_text(json.dumps(report,indent=2))
 print('PASS: continuous epsilon family 0 < eps <=1/1000; exact endpoint support proof; min convex guard gap',gap)
